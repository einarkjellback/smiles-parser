mod molecule;
mod molecule_building;

use molecule::*;
use molecule_building::*;
use std::iter::Peekable;

const BOND_TYPES: [(char, molecule::Order); 4] = [
    ('-', Order::Single),
    ('=', Order::Double),
    ('#', Order::Triple),
    ('$', Order::Quadruple),
];

const BRACKET_ATOM_START: char = '[';
const BRACKET_ATOM_END: char = ']';

const BRANCH_START: char = '(';
const BRANCH_END: char = ')';

pub fn build_from_smiles(smiles: &str) -> Result<FlatMol, String> {
    Ok(build_from_smiles_token(&mut smiles.chars().peekable())?.get_molecule_flatmol())
}

fn build_from_smiles_token<I>(tokens: &mut Peekable<I>) -> Result<SmilesFlatMol, String>
where
    I: Iterator<Item = char>,
{
    let mut mol_builder = SmilesFlatMol::new();
    let mut current_atom: Option<usize> = None; // Keeps track of which atoms to connect in the presence of branches
    let mut next_bond: Option<Bond> = None;
    while tokens.peek().is_some() {
        let token = tokens.next().unwrap();
        if token == BRACKET_ATOM_START {
            current_atom = if let Some(n) = current_atom {
                Some(n + 1)
            } else {
                Some(0)
            };

            if let Some(bond) = next_bond {
                mol_builder.add_bond(Bond {
                    order: bond.order,
                    from: bond.from,
                    to: current_atom.unwrap(),
                });
            }

            let hcount = add_bracket_atom(tokens, &mut mol_builder).unwrap();
            add_hydrogens(hcount, &mut mol_builder);

            let missing_bracket = String::from(format!(
                "Syntax error: Expected closing bracket '{}'",
                BRACKET_ATOM_END
            ));

            if let Some(c) = tokens.next() {
                if c != BRACKET_ATOM_END {
                    return Err(missing_bracket);
                }
            } else {
                return Err(missing_bracket);
            }

            let order = get_next_bond_order(tokens);
            next_bond = Some(Bond {
                order,
                from: current_atom.unwrap(),
                to: 0,
            });

            current_atom = if let Some(n) = current_atom {
                Some(n + hcount)
            } else {
                panic!("current_atom should not be None");
            };
        } else if token.is_ascii_uppercase() {
            let mut element = String::new();
            element.push(token);
            if let Some(c) = tokens.peek() {
                if c.is_ascii_lowercase() {
                    element.push(tokens.next().unwrap())
                }
            }
            let elem = Element::get_element(&element).unwrap();
            mol_builder.add_atom(SmilesAtom::new(
                Atom {
                    isotope: Isotope {
                        elem,
                        mass_number: MassNumber::Average,
                    },
                    charge: 0,
                    chiral: None,
                    class: None,
                },
                false,
            ));
        } else if token == BRANCH_START {
            let first_bond = tokens.peek().unwrap();
            let order = if let Some(&(_, order)) = BOND_TYPES.iter().find(|&(c, _)| c == first_bond)
            {
                tokens.next(); // Consume character denoting bond order
                order
            } else {
                Order::Single
            };
            mol_builder.add_bond(Bond {
                order,
                from: current_atom.unwrap(),
                to: mol_builder.atom_cnt(),
            });
            let branch = build_from_smiles_token(tokens)?;
            let offset = mol_builder.atom_cnt();
            for &atom in branch.atoms() {
                mol_builder.add_atom(atom);
            }
            for bond in branch.bonds() {
                mol_builder.add_bond(Bond {
                    order: bond.order,
                    from: bond.from + offset,
                    to: bond.to + offset,
                });
            }
        } else if token == BRANCH_END {
            return Ok(mol_builder);
        }
    }
    Ok(mol_builder)
}

fn add_bracket_atom<I>(
    tokens: &mut Peekable<I>,
    mol_builder: &mut SmilesFlatMol,
) -> Result<usize, String>
where
    I: Iterator<Item = char>,
{
    let mass_number = if let Ok(n) = read_number(tokens) {
        MassNumber::Exact(n)
    } else {
        MassNumber::Average
    };

    let mut elem = String::new();
    elem.push(tokens.next().unwrap());
    if tokens.peek().unwrap().is_ascii_lowercase() {
        elem.push(tokens.next().unwrap());
    }
    let elem = Element::get_element(&elem).unwrap();

    let chiral = if tokens.peek().unwrap() == &'@' {
        tokens.next(); // Consume '@' character
        Some(if tokens.peek().unwrap() == &'@' {
            tokens.next(); // Consume '@' character
            Chirality::Clockwise
        } else {
            Chirality::AntiClockwise
        })
    } else {
        None
    };

    let hcount = if tokens.peek().unwrap() == &'H' {
        tokens.next(); // Consume 'H' character
        if let Ok(n) = read_number(tokens) {
            n
        } else {
            1
        }
    } else {
        0
    };

    let c = tokens.peek().unwrap();
    let charge = if c == &'+' || c == &'-' {
        let sign: isize = if tokens.next().unwrap() == '+' { 1 } else { -1 };
        let magnitude = if let Ok(n) = read_number(tokens) {
            n as isize
        } else {
            1
        };
        sign * magnitude
    } else {
        0
    };

    let class = if tokens.peek().unwrap() == &':' {
        tokens.next(); // Consume ':' character
        Some(read_number(tokens)?)
    } else {
        None
    };

    mol_builder.add_atom(SmilesAtom::new(
        Atom {
            isotope: Isotope { elem, mass_number },
            charge,
            chiral,
            class,
        },
        false,
    ));
    Ok(hcount)
}

fn read_number<I>(tokens: &mut Peekable<I>) -> Result<usize, String>
where
    I: Iterator<Item = char>,
{
    if tokens.peek().unwrap().is_digit(10) {
        Ok((0..)
            .into_iter()
            .map(|_| {
                let c = tokens.peek();
                if c.unwrap().is_digit(10) {
                    Some(tokens.next().unwrap())
                } else {
                    None
                }
            })
            .take_while(|e| e != &None)
            .map(|e| e.unwrap().to_digit(10).unwrap())
            .collect::<Vec<_>>() // Collect so that we can reverse later
            .into_iter()
            .rev()
            .zip(0..)
            .map(|(digit, places)| digit as usize * usize::pow(10, places))
            .fold(0, |acc, e| acc + e))
    } else {
        Err(String::from("No digits detected."))
    }
}

fn add_hydrogens(hcount: usize, mol_builder: &mut SmilesFlatMol) {
    let from = mol_builder.atom_cnt() - 1;
    for _ in 0..hcount {
        let hydrogen = SmilesAtom::new(
            Atom {
                isotope: Isotope {
                    elem: Element::H,
                    mass_number: MassNumber::Average,
                },
                charge: 0,
                chiral: None,
                class: None,
            },
            false,
        );
        mol_builder.add_atom(hydrogen);
        mol_builder.add_bond(Bond {
            order: Order::Single,
            from,
            to: mol_builder.atom_cnt(),
        });
    }
}

fn get_next_bond_order<I>(tokens: &mut Peekable<I>) -> Order
where
    I: Iterator<Item = char>,
{
    let bond = tokens.peek();
    if let Some(bond) = bond {
        if let Some(&(_, order)) = BOND_TYPES.iter().find(|&&(c, _)| bond == &c) {
            tokens.next(); // Consume bond order character
            order
        } else {
            Order::Single
        }
    } else {
        Order::Single
    }
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn smiles_conv_single_letter_element() {
        let smiles = vec!["[C]", "[H]", "[B]"];
        let expected = vec![Element::C, Element::H, Element::B];
        for (smiles, expected) in smiles.into_iter().zip(expected.into_iter()) {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .isotope
                .elem;
            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_double_letter_element() {
        let smiles = vec!["[Co]", "[He]", "[Br]"];
        let expected = vec![Element::Co, Element::He, Element::Br];
        for (smiles, expected) in smiles.into_iter().zip(expected.into_iter()) {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .isotope
                .elem;
            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_diatomic_bonds() {
        let cases = vec![
            ("[C]-[C]", Order::Single),
            ("[C]=[C]", Order::Double),
            ("[C]#[C]", Order::Triple),
            ("[C]$[C]", Order::Quadruple),
        ];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .bonds[0]
                .order;

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_implicit_single_bond() {
        let smiles = "[C][C]";
        let expected = Order::Single;

        let actual = build_from_smiles(smiles)
            .expect(&format!("Testcase {} failed", smiles))
            .bonds[0]
            .order;

        assert_eq!(expected, actual, "Failed for case: {}", smiles);
    }

    #[test]
    fn smiles_conv_isotope() {
        let smiles_list = vec!["[14C]", "[014C]"];
        let expected = MassNumber::Exact(14);
        for smiles in smiles_list {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .isotope
                .mass_number;

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_charge() {
        let cases = vec![("[C+]", 1), ("[C+1]", 1), ("[C-]", -1), ("[C-15]", -15)];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .charge;

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_class() {
        let cases = vec![("[C:5]", 5), ("[C:005]", 5), ("[C:00]", 0)];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .class
                .expect(&format!("Testcase {} failed", smiles));

            assert_eq!(expected, actual);
        }
    }

    #[test]
    fn smiles_conv_chiral() {
        let cases = vec![
            ("[C@]", Chirality::AntiClockwise),
            ("[C@@]", Chirality::Clockwise),
        ];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0]
                .chiral
                .expect(&format!("Testcase {} failed", smiles));

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_hydrogens() {
        let cases = vec![
            ("[C]", 0),
            ("[CH0]", 0),
            ("[CH]", 1),
            ("[CH1]", 1),
            ("[CH9]", 9),
        ];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms
                .into_iter()
                .filter(|a| a.isotope.elem == Element::H)
                .count();

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_hydrogen_single_bond() {
        let smiles = "[CH]";
        let expected = Order::Single;

        let actual = build_from_smiles(smiles)
            .expect(&format!("Testcase {} failed", smiles))
            .bonds[0]
            .order;

        assert_eq!(expected, actual);
    }

    #[test]
    fn smiles_conv_hydrogen_from() {
        let smiles = "[CH4]";
        let expected = 0;

        for bond in build_from_smiles(smiles)
            .expect(&format!("Testcase {} failed", smiles))
            .bonds
        {
            let actual = bond.from;
            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_organic_subset_atoms() {
        use Element::*;

        let mass_number = MassNumber::Average;
        let charge = 0;
        let chiral = None;
        let class = None;

        let cases = vec![
            ("B", B),
            ("C", C),
            ("N", N),
            ("O", O),
            ("P", P),
            ("S", S),
            ("F", F),
            ("Cl", Cl),
            ("Br", Br),
            ("I", I),
        ]
        .into_iter()
        .map(|(s, elem)| (s, Isotope { elem, mass_number }))
        .map(|(s, isotope)| {
            (
                s,
                Atom {
                    isotope,
                    charge,
                    chiral,
                    class,
                },
            )
        })
        .collect::<Vec<_>>();
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms[0];

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[test]
    fn smiles_conv_organic_subset_implicit_hydrogens() {
        let cases = vec![
            ("B", 3),
            ("C", 4),
            ("N", 3),
            ("[H]N([H])([H])[H]", 5),
            ("O", 2),
            ("P", 3),
            ("[H]P([H])([H])[H]", 5),
            ("S", 2),
            ("[H]S([H])[H]", 4),
            ("[H]S([H])([H])([H])[H]", 6),
            ("F", 1),
            ("Cl", 1),
            ("Br", 1),
            ("I", 1),
        ];
        for (smiles, expected) in cases {
            let actual = build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .atoms
                .into_iter()
                .filter(|a| a.isotope.elem == Element::H)
                .count();

            assert_eq!(expected, actual, "Failed for case: {}", smiles);
        }
    }

    #[ignore]
    #[test]
    fn smiles_conv_branch_neighbors() {
        let cases = vec![
            "[C](-[H])",
            "[C]-(-[H])[H]",
            "[C](-[H])(-[H])",
            "[C]-(-[H])(-[H])[H]",
        ];
        for smiles in cases {
            assert!(build_from_smiles(smiles)
                .expect(&format!("Testcase {} failed", smiles))
                .bonds
                .iter()
                .all(|b| b.from == 0));
        }
    }
}
