use crate::molecule::*;

#[derive(Clone, Copy)]
pub struct SmilesAtom {
    atom: Atom,
    from_organic_subset: bool,
}

impl SmilesAtom {
    pub fn new(atom: Atom, from_organic_subset: bool) -> SmilesAtom {
        SmilesAtom {
            atom,
            from_organic_subset,
        }
    }

    pub fn get_atom(self) -> Atom {
        self.atom
    }

    fn get_missing_hydrogens(&self, neighbor_count: usize) -> usize {
        use Element::*;
        if self.from_organic_subset {
            let n = neighbor_count;
            match self.atom.isotope.elem {
                F | Cl | Br | I => 1,
                O => 2,
                B => 3,
                C => 4,
                N | P => {
                    if n <= 3 {
                        3 - n
                    } else {
                        usize::max(0, 5 - n)
                    }
                }
                S => {
                    if n <= 2 {
                        2 - n
                    } else if n <= 4 {
                        4 - n
                    } else {
                        usize::max(0, 6 - n)
                    }
                }
                _ => 0,
            }
        } else {
            0
        }
    }
}

pub struct SmilesFlatMol {
    atoms: Vec<SmilesAtom>,
    bonds: Vec<Bond>,
}

impl SmilesFlatMol {
    pub fn new() -> SmilesFlatMol {
        SmilesFlatMol {
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    pub fn add_atom(&mut self, atom: SmilesAtom) {
        self.atoms.push(atom);
    }

    pub fn atom_cnt(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &Vec<SmilesAtom> {
        &self.atoms
    }

    pub fn bonds(&self) -> &Vec<Bond> {
        &self.bonds
    }

    pub fn get_molecule_flatmol(self) -> FlatMol {
        FlatMol {
            atoms: self.atoms.into_iter().map(|a| a.get_atom()).collect(),
            bonds: self.bonds,
        }
    }

    pub fn add_bond(&mut self, bond: Bond) {
        let i = self.add_bond_at(bond, 0, self.bonds.len());
        self.bonds.insert(i, bond);
    }

    fn add_bond_at(&self, bond: Bond, lo: usize, hi: usize) -> usize {
        if hi == lo {
            lo
        } else {
            let mid = lo + (hi - lo) / 2;
            let other = self.bonds[mid];
            if other.from < bond.from {
                self.add_bond_at(bond, mid + 1, hi)
            } else if other.from > bond.from {
                self.add_bond_at(bond, lo, mid)
            } else if other.to < bond.to {
                self.add_bond_at(bond, mid + 1, hi)
            } else {
                self.add_bond_at(bond, lo, mid)
            }
        }
    }
}
