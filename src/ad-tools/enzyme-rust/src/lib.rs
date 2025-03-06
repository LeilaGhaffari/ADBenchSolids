#![allow(non_snake_case)]
#![cfg_attr(feature = "enzyme", feature(autodiff))]

use std::ops::{Add, Mul, Sub};

type Mat3x3 = [[f64; 3]; 3];

fn pack_mat(input: [f64; 9]) -> Mat3x3 {
    let mut out: Mat3x3 = [[0.; 3]; 3];
    out[0][0] = input[0];
    out[0][1] = input[1];
    out[0][2] = input[2];
    out[1][0] = input[3];
    out[1][1] = input[4];
    out[1][2] = input[5];
    out[2][0] = input[6];
    out[2][1] = input[7];
    out[2][2] = input[8];
    out
}

fn matmul(a: &Mat3x3, at: bool, b: &Mat3x3, bt: bool) -> Mat3x3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let aa = if at { a[k][i] } else { a[i][k] };
                let bb = if bt { b[j][k] } else { b[k][j] };
                c[i][j] += aa * bb;
            }
        }
    }
    c
}

fn matadd(alpha: f64, a: &Mat3x3, beta: f64, b: &Mat3x3) -> Mat3x3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = alpha * a[i][j] + beta * b[i][j];
        }
    }
    c
}

fn matinv(a: &Mat3x3) -> Mat3x3 {
    let mut b = [
        [
            a[1][1] * a[2][2] - a[1][2] * a[2][1],
            a[0][2] * a[2][1] - a[0][1] * a[2][2],
            a[0][1] * a[1][2] - a[0][2] * a[1][1],
        ],
        [
            a[1][2] * a[2][0] - a[1][0] * a[2][2],
            a[0][0] * a[2][2] - a[0][2] * a[2][0],
            a[0][2] * a[1][0] - a[0][0] * a[1][2],
        ],
        [
            a[1][0] * a[2][1] - a[1][1] * a[2][0],
            a[0][1] * a[2][0] - a[0][0] * a[2][1],
            a[0][0] * a[1][1] - a[0][1] * a[1][0],
        ],
    ];
    let det = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] /= det;
        }
    }
    b
}

fn deformation_gradient(grad_u: Mat3x3) -> Mat3x3 {
    let mut F = grad_u;
    for i in 0..3 {
        F[i][i] = F[i][i] + 1.;
    }
    F
}

/// Kelvin-Mandel 6-vector representation of symmetric 3x3 matrices.
///
/// [Kelvin-Mandel](https://en.wikipedia.org/wiki/Voigt_notation#Mandel_notation)
/// notation has the some ordering as Voigt, but a normalization such that the
/// Frobenius matrix inner product is equivalent to the standard inner product
/// on KM vectors. This is useful for geometry in return mapping algorithms for
/// plasticity.
#[derive(Debug, Clone, Copy)]
pub struct KM {
    pub vals: [f64; 6],
}

impl KM {
    pub fn from_voigt(v: &[f64; 6]) -> Self {
        let d = 2.0_f64.sqrt();
        Self {
            vals: [v[0], v[1], v[2], d * v[3], d * v[4], d * v[5]],
        }
    }

    pub fn zero() -> Self {
        Self { vals: [0.0; 6] }
    }

    pub fn identity() -> Self {
        Self {
            vals: [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        }
    }

    pub fn scale(&self, scale: f64) -> Self {
        Self {
            vals: self.vals.map(|v| v * scale),
        }
    }

    pub fn from_binary_op(x: KM, y: KM, op: impl Fn(f64, f64) -> f64) -> Self {
        let mut result = KM::zero();
        result
            .vals
            .iter_mut()
            .zip(x.vals.into_iter().zip(y.vals.into_iter()))
            .for_each(|(r, (x, y))| *r = op(x, y));
        result
    }

    pub fn from_matrix(a: Mat3x3) -> Self {
        let s = 0.5_f64.sqrt();
        KM {
            vals: [
                a[0][0],
                a[1][1],
                a[2][2],
                s * (a[1][2] + a[2][1]),
                s * (a[0][2] + a[2][0]),
                s * (a[0][1] + a[1][0]),
            ],
        }
    }

    pub fn to_matrix(&self) -> Mat3x3 {
        let d = 0.5_f64.sqrt();
        let v = &self.vals;
        [
            [v[0], d * v[5], d * v[4]],
            [d * v[5], v[1], d * v[3]],
            [d * v[4], d * v[3], v[2]],
        ]
    }

    pub fn to_voigt(&self) -> [f64; 6] {
        let d = 2.0_f64.sqrt();
        let v = &self.vals;
        [v[0], v[1], v[2], v[3] / d, v[4] / d, v[5] / d]
    }

    // Stably evaluate e = (F F^T - I)/2 from displacement gradient H = F - I.
    pub fn green_euler(H: Mat3x3) -> Self {
        Self::from_matrix(H) + 0.5 * Self::from_matrix(matmul(&H, false, &H, true))
    }

    // infinitesimal strain using current configuration displacement gradient ddu/dx = ddu/dX F^{-1}
    pub fn epsilon(ddudx: &Mat3x3) -> Self {
        Self::from_matrix(*ddudx)
    }

    pub fn trace(&self) -> f64 {
        let v = &self.vals;
        v[0] + v[1] + v[2]
    }

    pub fn det(&self) -> f64 {
        let m = self.to_matrix();
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    // Compute b = I + 2e or C = I + 2E.
    pub fn cauchy_green(&self) -> Self {
        let v = &self.vals;
        Self {
            vals: [
                1.0 + 2.0 * v[0],
                1.0 + 2.0 * v[1],
                1.0 + 2.0 * v[2],
                2.0 * v[3],
                2.0 * v[4],
                2.0 * v[5],
            ],
        }
    }

    pub fn inv(&self) -> Self {
        let a = self.to_matrix();
        KM::from_matrix(matinv(&a))
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.vals
            .into_iter()
            .zip(other.vals)
            .map(|(a, b)| a * b)
            .sum()
    }

    pub fn norm(&self) -> f64 {
        self.dot(self).sqrt()
    }
}

impl Add for KM {
    type Output = KM;
    fn add(self, rhs: Self) -> Self::Output {
        KM::from_binary_op(self, rhs, f64::add)
    }
}

impl Mul for KM {
    type Output = KM;
    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.to_matrix();
        let b = rhs.to_matrix();
        let c = matmul(&a, false, &b, false);
        KM::from_matrix(c)
    }
}

impl Mul<&KM> for f64 {
    type Output = KM;
    fn mul(self, rhs: &KM) -> Self::Output {
        rhs.scale(self)
    }
}

impl Mul<KM> for f64 {
    type Output = KM;
    fn mul(self, rhs: KM) -> Self::Output {
        rhs.scale(self)
    }
}

impl Sub for KM {
    type Output = KM;
    fn sub(self, rhs: Self) -> Self::Output {
        KM::from_binary_op(self, rhs, f64::sub)
    }
}

pub struct NH {
    pub mu: f64,
    pub lambda: f64,
}

impl NH {
    pub fn from_lame(lambda: f64, mu: f64) -> Self {
        Self { lambda, mu }
    }
    // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    pub fn from_youngs(E: f64, nu: f64) -> Self {
        Self::from_lame(
            E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
            E / (2.0 * (1.0 + nu)),
        )
    }

    #[cfg(feature = "enzyme")]
    pub fn stress(&self, e: &KM) -> KM {
        let mut tau = KM::zero();
        enzyme::stress_enz(e, self, &mut tau);
        tau
    }

    #[cfg(feature = "enzyme")]
    pub fn d_stress(&self, e: &KM, de: &KM) -> KM {
        let mut tau = KM::zero();
        let mut dtau = KM::zero();
        enzyme::d_stress_enz(e, de, self, &mut tau, &mut dtau);
        dtau
    }
}

// We can only differentiate free functions, not methods (yet)
// Helmholtz free energy density
#[cfg(feature = "enzyme")]
mod enzyme {
    use crate::{KM, NH};
    use std::autodiff::autodiff;

    #[autodiff(d_psi, Reverse, Duplicated, Const, Active)]
    pub fn psi(e: &KM, nh: &NH) -> f64 {
        let mu = nh.mu;
        let lambda = nh.lambda;
        let J = e.cauchy_green().det().sqrt();
        0.25 * lambda * (J * J - 1.0 - 2.0 * J.ln()) + mu * (e.trace() - J.ln())
    }

    #[autodiff(d_stress_enz, Forward, Dual, Const, Dual)]
    pub fn stress_enz(e: &KM, nh: &NH, tau: &mut KM) {
        let mut dpsi_de = KM::zero();
        d_psi(&e, &mut dpsi_de, &nh, 1.0);
        let b = e.cauchy_green();
        *tau = dpsi_de * b;
    }
}

pub mod analytic {
    use super::*;
    pub fn stress(e: &KM, nh: &NH) -> KM {
        let lambda = nh.lambda;
        let mu = nh.mu;
        let b = e.cauchy_green();
        let J = b.det().sqrt();
        let I = KM::identity();
        0.5 * lambda * (J * J - 1.0) * I + 2.0 * mu * e
    }

    pub fn d_stress(e: &KM, ddudx: &Mat3x3, nh: &NH) -> KM {
        let lambda = nh.lambda;
        let mu = nh.mu;
        let b = e.cauchy_green();
        let J = b.det().sqrt();
        let tau = stress(e, nh).to_matrix();
        let deps = KM::epsilon(ddudx);
        let I = KM::identity();
        let FdSFt = lambda * J * J * deps.trace() * I + (2.0 * mu - lambda * (J * J - 1.0)) * deps;
        2.0 * KM::from_matrix(matmul(&ddudx, false, &tau, false)) + FdSFt
    }
}

fn q_data_pack_mat(Q: usize, i: usize, input: &[f64]) -> Mat3x3 {
    let mut out: Mat3x3 = [[0.; 3]; 3];
    out[0][0] = input[0 * Q + i];
    out[0][1] = input[1 * Q + i];
    out[0][2] = input[2 * Q + i];
    out[1][0] = input[3 * Q + i];
    out[1][1] = input[4 * Q + i];
    out[1][2] = input[5 * Q + i];
    out[2][0] = input[6 * Q + i];
    out[2][1] = input[7 * Q + i];
    out[2][2] = input[8 * Q + i];
    out
}

fn q_data_unpack_mat(Q: usize, i: usize, input: Mat3x3, out: &mut [f64]) {
    out[0 * Q + i] = input[0][0];
    out[1 * Q + i] = input[0][1];
    out[2 * Q + i] = input[0][2];
    out[3 * Q + i] = input[1][0];
    out[4 * Q + i] = input[1][1];
    out[5 * Q + i] = input[1][2];
    out[6 * Q + i] = input[2][0];
    out[7 * Q + i] = input[2][1];
    out[8 * Q + i] = input[2][2];
}

fn stored_values_pack(
    Q: usize,
    i: usize,
    start: usize,
    num_comp: usize,
    local: &[f64],
    stored: &mut [f64],
) {
    for j in 0..num_comp {
        stored[(start + j) * Q + i] = local[j];
    }
}

#[allow(dead_code)]
fn stored_values_unpack(
    Q: usize,
    i: usize,
    start: usize,
    num_comp: usize,
    stored: &[f64],
    local: &mut [f64],
) {
    for j in 0..num_comp {
        local[j] = stored[(start + j) * Q + i];
    }
}

#[no_mangle]
pub extern "C" fn compute_f_analytic(
    Q: usize,
    mu: f64,
    lambda: f64,
    dXdx_init: *const f64,
    dudX: *const f64,
    num_stored_comp: usize,
    stored_values: *mut f64,
    f1: *mut f64,
) {
    let dXdx_init = unsafe { std::slice::from_raw_parts(dXdx_init, 9 * Q) };
    let dudX = unsafe { std::slice::from_raw_parts(dudX, 9 * Q) };
    let stored_values =
        unsafe { std::slice::from_raw_parts_mut(stored_values, num_stored_comp * Q) };
    let f1 = unsafe { std::slice::from_raw_parts_mut(f1, 9 * Q) };

    let nh = NH::from_lame(lambda, mu);
    for i in 0..Q {
        let dXdx_init_loc = q_data_pack_mat(Q, i, dXdx_init);
        let dudX_loc = q_data_pack_mat(Q, i, dudX);
        let Grad_u = matmul(&dudX_loc, false, &dXdx_init_loc, false);
        let F = deformation_gradient(Grad_u);
        let Finv = matinv(&F);
        let dXdx = matmul(&dXdx_init_loc, false, &Finv, false);
        let e_sym = KM::green_euler(Grad_u);
        let tau_sym = analytic::stress(&e_sym, &nh);
        q_data_unpack_mat(Q, i, tau_sym.to_matrix(), f1);

        let dXdx_flat: Vec<f64> = dXdx.iter().flatten().copied().collect();
        stored_values_pack(Q, i, 0, 9, &dXdx_flat, stored_values);
        stored_values_pack(Q, i, 9, 6, &e_sym.vals, stored_values);
    }
}

#[no_mangle]
pub extern "C" fn compute_df_analytic(
    Q: usize,
    mu: f64,
    lambda: f64,
    ddudX: *const f64,
    num_stored_comp: usize,
    stored_values: *mut f64,
    df: *mut f64,
) {
    let ddudX = unsafe { std::slice::from_raw_parts(ddudX, 9 * Q) };
    let stored_values =
        unsafe { std::slice::from_raw_parts_mut(stored_values, num_stored_comp * Q) };
    let df = unsafe { std::slice::from_raw_parts_mut(df, 9 * Q) };

    let nh = NH::from_lame(lambda, mu);
    for i in 0..Q {
        let mut e_sym = KM::zero();
        let ddudX_loc = q_data_pack_mat(Q, i, ddudX);
        let mut dXdx_flat: [f64; 9] = [0.0; 9];
        stored_values_unpack(Q, i, 0, 9, stored_values, &mut dXdx_flat);
        stored_values_unpack(Q, i, 9, 6, stored_values, &mut e_sym.vals);
        let dXdx = pack_mat(dXdx_flat);
        let ddudx = matmul(&ddudX_loc, false, &dXdx, false);
        let tau_sym = analytic::stress(&e_sym, &nh);
        let tau = tau_sym.to_matrix();
        let ddudx_tau = matmul(&ddudx, false, &tau, false);
        let b = e_sym.cauchy_green();
        let J = b.det().sqrt();
        let deps = KM::epsilon(&ddudx);
        let I = KM::identity();
        let FdSFt = lambda * J * J * deps.trace() * I + (2.0 * mu - lambda * (J * J - 1.0)) * deps;
        let df_mat = matadd(1., &ddudx_tau, 1., &FdSFt.to_matrix());
        q_data_unpack_mat(Q, i, df_mat, df);
    }
}

#[no_mangle]
#[cfg(feature = "enzyme")]
pub extern "C" fn compute_f_enzyme(
    Q: usize,
    mu: f64,
    lambda: f64,
    dXdx_init: *const f64,
    dudX: *const f64,
    num_stored_comp: usize,
    stored_values: *mut f64,
    f1: *mut f64,
) {
    let dXdx_init = unsafe { std::slice::from_raw_parts(dXdx_init, 9 * Q) };
    let dudX = unsafe { std::slice::from_raw_parts(dudX, 9 * Q) };
    let stored_values =
        unsafe { std::slice::from_raw_parts_mut(stored_values, num_stored_comp * Q) };
    let f1 = unsafe { std::slice::from_raw_parts_mut(f1, 9 * Q) };

    let nh = NH::from_lame(lambda, mu);
    for i in 0..Q {
        let dXdx_init_loc = q_data_pack_mat(Q, i, dXdx_init);
        let dudX_loc = q_data_pack_mat(Q, i, dudX);
        let Grad_u = matmul(&dudX_loc, false, &dXdx_init_loc, false);
        let F = deformation_gradient(Grad_u);
        let Finv = matinv(&F);
        let dXdx = matmul(&dXdx_init_loc, false, &Finv, false);
        let e_sym = KM::green_euler(Grad_u);
        let mut tau_sym = KM::zero();
        enzyme::stress_enz(&e_sym, &nh, &mut tau_sym);
        q_data_unpack_mat(Q, i, tau_sym.to_matrix(), f1);

        let dXdx_flat: Vec<f64> = dXdx.iter().flatten().copied().collect();
        stored_values_pack(Q, i, 0, 9, &dXdx_flat, stored_values);
        stored_values_pack(Q, i, 9, 6, &e_sym.vals, stored_values);
    }
}

#[no_mangle]
#[cfg(feature = "enzyme")]
pub extern "C" fn compute_df_enzyme(
    Q: usize,
    mu: f64,
    lambda: f64,
    ddudX: *const f64,
    num_stored_comp: usize,
    stored_values: *mut f64,
    df: *mut f64,
) {
    let ddudX = unsafe { std::slice::from_raw_parts(ddudX, 9 * Q) };
    let stored_values =
        unsafe { std::slice::from_raw_parts_mut(stored_values, num_stored_comp * Q) };
    let df = unsafe { std::slice::from_raw_parts_mut(df, 9 * Q) };

    let nh = NH::from_lame(lambda, mu);
    for i in 0..Q {
        let mut e_sym = KM::zero();
        let ddudX_loc = q_data_pack_mat(Q, i, ddudX);
        let mut dXdx_flat: [f64; 9] = [0.0; 9];
        stored_values_unpack(Q, i, 0, 9, stored_values, &mut dXdx_flat);
        stored_values_unpack(Q, i, 9, 6, stored_values, &mut e_sym.vals);
        let dXdx = pack_mat(dXdx_flat);
        let ddudx = matmul(&ddudX_loc, false, &dXdx, false);
        let tau_sym = analytic::stress(&e_sym, &nh);
        let tau = tau_sym.to_matrix();
        let deps = KM::epsilon(&ddudx);
        let de_sym = 2.0 * KM::from_matrix(matmul(&ddudx, false, &e_sym.to_matrix(), false)) + deps;
        let dtau_sym = nh.d_stress(&e_sym, &de_sym);
        let dtau = dtau_sym.to_matrix();
        let tau_ddudx = matmul(&tau, false, &ddudx, true);
        let df_mat = matadd(1., &dtau, -1., &tau_ddudx);
        q_data_unpack_mat(Q, i, df_mat, df);
    }
}
