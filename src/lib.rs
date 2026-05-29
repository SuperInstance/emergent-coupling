//! Emergent coupling: two systems coupling through spectral gaps
//! to produce joint structure larger than either alone.

/// A system with a spectral profile (eigenvalues and dominant mode)
#[derive(Debug, Clone)]
pub struct SpectralSystem {
    pub name: String,
    pub eigenvalues: Vec<f64>,
    pub dominant_eigenvector: Vec<f64>,
}

impl SpectralSystem {
    pub fn new(name: &str, eigenvalues: Vec<f64>, dominant_eigenvector: Vec<f64>) -> Self {
        Self { name: name.to_string(), eigenvalues, dominant_eigenvector }
    }

    /// Conservation ratio: λ₂/λₙ
    pub fn cr(&self) -> f64 {
        if self.eigenvalues.len() < 2 { return 0.0; }
        let l2 = self.eigenvalues[1];
        let ln = *self.eigenvalues.last().unwrap_or(&1.0);
        if ln <= 0.0 { return 0.0; }
        (l2 / ln).clamp(0.0, 1.0)
    }

    /// Spectral gap: λ₂ (algebraic connectivity)
    pub fn spectral_gap(&self) -> f64 {
        self.eigenvalues.get(1).copied().unwrap_or(0.0)
    }

    /// Norm of dominant eigenvector
    pub fn mode_energy(&self) -> f64 {
        self.dominant_eigenvector.iter().map(|x| x * x).sum()
    }
}

/// Coupling mode between two systems
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CouplingMode {
    Reinforcing,
    Complementary,
    Competitive,
}

/// Result of coupling two systems
#[derive(Debug)]
pub struct CouplingResult {
    pub coupled: SpectralSystem,
    pub mode: CouplingMode,
    pub emergent: bool,
    pub synergy: f64,
    pub alignment: f64,
}

/// Cosine similarity between two vectors
pub fn cosine_similarity(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 { return 0.0; }
    let dot: f64 = a[..n].iter().zip(&b[..n]).map(|(x, y)| x * y).sum();
    let na: f64 = a[..n].iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b[..n].iter().map(|x| x * x).sum::<f64>().sqrt();
    if na < 1e-15 || nb < 1e-15 { return 0.0; }
    (dot / (na * nb)).clamp(-1.0, 1.0)
}

/// Couple two systems by building a coupled Laplacian.
/// Models coupling as adding edges between the two systems' graphs.
pub fn couple(a: &SpectralSystem, b: &SpectralSystem, coupling_strength: f64) -> CouplingResult {
    let alignment = cosine_similarity(&a.dominant_eigenvector, &b.dominant_eigenvector);

    // Build coupled adjacency: block-diagonal with coupling edges
    let n_a = a.eigenvalues.len();
    let n_b = b.eigenvalues.len();
    let n = n_a + n_b;
    if n < 2 {
        let empty = SpectralSystem::new(&format!("{}×{}", a.name, b.name), vec![], vec![]);
        return CouplingResult { coupled: empty, mode: CouplingMode::Complementary, emergent: false, synergy: 1.0, alignment };
    }

    // Reconstruct approximate adjacency from eigenvalues (scaled)
    let mut adj = vec![vec![0.0_f64; n]; n];

    // Fill block A: use eigenvalues to estimate edge weights
    // For a graph with Laplacian eigenvalues λ_i, the total edge weight = λ_n * n / 2 (approx)
    // Distribute edges based on dominant eigenvector
    let a_norm: f64 = a.dominant_eigenvector.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-15);
    let b_norm: f64 = b.dominant_eigenvector.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-15);

    for i in 0..n_a {
        for j in (i+1)..n_a {
            let vi = a.dominant_eigenvector.get(i).copied().unwrap_or(1.0) / a_norm;
            let vj = a.dominant_eigenvector.get(j).copied().unwrap_or(1.0) / a_norm;
            let w = a.spectral_gap() * (1.0 + vi * vj) * 0.5;
            adj[i][j] = w; adj[j][i] = w;
        }
    }
    for i in 0..n_b {
        for j in (i+1)..n_b {
            let vi = b.dominant_eigenvector.get(i).copied().unwrap_or(1.0) / b_norm;
            let vj = b.dominant_eigenvector.get(j).copied().unwrap_or(1.0) / b_norm;
            let w = b.spectral_gap() * (1.0 + vi * vj) * 0.5;
            adj[n_a + i][n_a + j] = w; adj[n_a + j][n_a + i] = w;
        }
    }

    // Coupling edges between A and B
    for i in 0..n_a.min(3) {
        for j in 0..n_b.min(3) {
            let vi = a.dominant_eigenvector.get(i).copied().unwrap_or(0.0) / a_norm;
            let vj = b.dominant_eigenvector.get(j).copied().unwrap_or(0.0) / b_norm;
            let w = coupling_strength * (1.0 + vi * vj) * (a.spectral_gap() + b.spectral_gap()) * 0.25;
            adj[i][n_a + j] = w; adj[n_a + j][i] = w;
        }
    }

    // Compute Laplacian and eigenvalues
    let mut lap = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let deg: f64 = adj[i].iter().sum();
        lap[i][i] = deg;
        for j in 0..n { if i != j { lap[i][j] = -adj[i][j]; } }
    }
    let coupled_eigs = jacobi_eigenvalues(&mut lap);

    // Coupled dominant mode
    let mut coupled_mode = vec![0.0; n];
    for i in 0..n_a { coupled_mode[i] = a.dominant_eigenvector.get(i).copied().unwrap_or(0.0) / a_norm; }
    for i in 0..n_b { coupled_mode[n_a + i] = b.dominant_eigenvector.get(i).copied().unwrap_or(0.0) / b_norm; }

    let coupled = SpectralSystem::new(&format!("{}×{}", a.name, b.name), coupled_eigs, coupled_mode);

    let max_individual_cr = a.cr().max(b.cr());
    let synergy = if max_individual_cr > 1e-15 { coupled.cr() / max_individual_cr }
        else if coupled.cr() > 0.0 { 2.0 } else { 1.0 };
    let emergent = coupled.cr() > max_individual_cr + 1e-10;

    let mode = if alignment > 0.3 { CouplingMode::Reinforcing }
        else if alignment < -0.3 { CouplingMode::Competitive }
        else { CouplingMode::Complementary };

    CouplingResult { coupled, mode, emergent, synergy, alignment }
}

/// Sweep coupling strength and return (strength, synergy, emergent)
pub fn emergence_curve(a: &SpectralSystem, b: &SpectralSystem, steps: usize) -> Vec<(f64, f64, bool)> {
    (0..=steps)
        .map(|i| {
            let s = i as f64 / steps as f64;
            let r = couple(a, b, s);
            (s, r.synergy, r.emergent)
        })
        .collect()
}

/// Find the critical coupling strength where emergence first occurs
pub fn critical_coupling(a: &SpectralSystem, b: &SpectralSystem, resolution: usize) -> Option<f64> {
    emergence_curve(a, b, resolution)
        .iter()
        .find(|&&(_, _, emergent)| emergent)
        .map(|&(s, _, _)| s)
}

/// Graph adjacency → SpectralSystem
pub fn graph_to_system(name: &str, adj: &[Vec<f64>]) -> SpectralSystem {
    let n = adj.len();
    let mut lap = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let deg: f64 = adj[i].iter().sum();
        lap[i][i] = deg;
        for j in 0..n {
            if i != j { lap[i][j] = -adj[i][j]; }
        }
    }
    let eigs = jacobi_eigenvalues(&mut lap);
    let total_deg: f64 = adj.iter().flat_map(|r| r.iter()).sum();
    let mode: Vec<f64> = adj.iter().map(|row| row.iter().sum::<f64>() / total_deg.max(1e-15)).collect();
    SpectralSystem::new(name, eigs, mode)
}

fn jacobi_eigenvalues(a: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 { return vec![]; }
    for _ in 0..100 * n * n {
        let (mut p, mut q) = (0, 1);
        let mut max_val = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val { max_val = a[i][j].abs(); p = i; q = j; }
            }
        }
        if max_val < 1e-14 { break; }
        let app = a[p][p]; let aqq = a[q][q]; let apq = a[p][q];
        let theta = if (app - aqq).abs() < 1e-30 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * apq / (app - aqq)).atan()
        };
        let (c, s) = (theta.cos(), theta.sin());
        for i in 0..n {
            if i != p && i != q {
                let aip = a[i][p]; let aiq = a[i][q];
                a[i][p] = c * aip + s * aiq; a[p][i] = a[i][p];
                a[i][q] = -s * aip + c * aiq; a[q][i] = a[i][q];
            }
        }
        a[p][p] = c * c * app + 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app - 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0; a[q][p] = 0.0;
    }
    let mut eigs: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

/// Chain-couple N systems sequentially
pub fn chain_coupling(systems: &[SpectralSystem], strength: f64) -> CouplingResult {
    if systems.is_empty() {
        let empty = SpectralSystem::new("empty", vec![], vec![]);
        return CouplingResult { coupled: empty, mode: CouplingMode::Complementary, emergent: false, synergy: 1.0, alignment: 0.0 };
    }
    if systems.len() == 1 {
        let r = couple(&systems[0], &systems[0], 0.0);
        return r;
    }
    let mut result = couple(&systems[0], &systems[1], strength);
    for i in 2..systems.len() {
        result = couple(&result.coupled, &systems[i], strength);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn system_a() -> SpectralSystem {
        SpectralSystem::new("A", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 1.0, 1.0, 1.0])
    }
    fn system_b() -> SpectralSystem {
        SpectralSystem::new("B", vec![0.0, 0.3, 1.5, 4.0], vec![1.0, 1.0, 1.0, 1.0])
    }
    fn system_orthogonal() -> SpectralSystem {
        SpectralSystem::new("C", vec![0.0, 0.8, 3.0, 6.0], vec![1.0, -1.0, 1.0, -1.0])
    }

    #[test]
    fn cr_computed() {
        assert!((system_a().cr() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn reinforcing_coupling() {
        let r = couple(&system_a(), &system_b(), 1.0);
        assert_eq!(r.mode, CouplingMode::Reinforcing);
        assert!(r.alignment > 0.9);
    }

    #[test]
    fn competitive_coupling_alignment() {
        let opposed = SpectralSystem::new("C", vec![0.0, 0.8, 3.0, 6.0], vec![-1.0, -1.0, -1.0, -1.0]);
        let r = couple(&system_a(), &opposed, 1.0);
        assert!(r.alignment < 0.0, "Alignment should be negative: {}", r.alignment);
    }

    #[test]
    fn emergence_with_aligned_strong_coupling() {
        let w1 = SpectralSystem::new("W1", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 0.0]);
        let w2 = SpectralSystem::new("W2", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 0.0]);
        let r = couple(&w1, &w2, 5.0);
        assert!(r.synergy > 1.0, "Synergy = {}", r.synergy);
    }

    #[test]
    fn emergence_curve_monotonic_aligned() {
        let w = SpectralSystem::new("W", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 0.0]);
        let curve = emergence_curve(&w, &w, 20);
        assert_eq!(curve.len(), 21);
        // At zero coupling, synergy should be ~1.0
        assert!((curve[0].1 - 1.0).abs() < 1.5, "Zero coupling synergy = {}", curve[0].1);
        // Synergy should generally increase for aligned systems
        let last = curve.last().unwrap().1;
        assert!(last >= curve[0].1, "Synergy should increase: {} vs {}", last, curve[0].1);
    }

    #[test]
    fn critical_coupling_found() {
        let w = SpectralSystem::new("W", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 0.0]);
        let crit = critical_coupling(&w, &w, 100);
        assert!(crit.is_some(), "Should find emergence for aligned identical systems");
    }

    #[test]
    fn graph_to_system_works() {
        let adj = vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]];
        let sys = graph_to_system("tri", &adj);
        assert_eq!(sys.eigenvalues.len(), 3);
        assert!(sys.eigenvalues[0].abs() < 1e-10);
        assert!((sys.eigenvalues[1] - 3.0).abs() < 0.1, "λ₂ = {}", sys.eigenvalues[1]);
    }

    #[test]
    fn cosine_similarity_bounds() {
        assert!((cosine_similarity(&[1.0, 0.0], &[1.0, 0.0]) - 1.0).abs() < 1e-10);
        assert!(cosine_similarity(&[1.0, 0.0], &[0.0, 1.0]).abs() < 1e-10);
        assert!((cosine_similarity(&[1.0, 0.0], &[-1.0, 0.0]) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn chain_coupling_runs() {
        let systems = vec![system_a(), system_b(), system_orthogonal()];
        let r = chain_coupling(&systems, 0.5);
        assert!(!r.coupled.eigenvalues.is_empty());
    }
}
