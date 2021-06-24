use std::{collections::HashMap, vec::Vec};
use Element::*;

pub struct FlatMol {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Atom {
    pub isotope: Isotope,
    pub charge: Charge,
    pub chiral: Option<Chirality>,
    pub class: Option<Class>,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Isotope {
    pub elem: Element,
    pub mass_number: MassNumber,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Element {
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Fl,
    Lv,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Ac,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
}

impl Element {
    pub fn get_element(element: &str) -> Result<Element, String> {
        let keys = vec![
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
            "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
            "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
            "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf", "Ta", "W", "Re",
            "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Rf",
            "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl", "Lv", "La", "Ce", "Pr", "Nd",
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa",
            "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        ];
        let elements = vec![
            H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr,
            Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd,
            Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi,
            Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Fl, Lv, La, Ce, Pr, Nd, Pm, Sm,
            Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md,
            No, Lr,
        ];

        let str_to_el: HashMap<_, _> = keys.into_iter().zip(elements.into_iter()).collect();
        match str_to_el.get(element) {
            Some(&e) => Ok(e),
            None => {
                let msg = format!("Element '{}' does not exist", element);
                Err(msg)
            }
        }
    }
}

type Charge = isize;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Chirality {
    Clockwise,
    AntiClockwise,
}

type Class = usize;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum MassNumber {
    Exact(usize),
    Average,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Bond {
    pub order: Order,
    pub from: usize,
    pub to: usize,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Order {
    Single,
    Double,
    Triple,
    Quadruple,
}
