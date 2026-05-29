//! Integration tests for emergent-coupling

use emergent_coupling::*;

fn simple_system(name: &str) -> SpectralSystem {
    SpectralSystem::new(name, vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 1.0, 1.0, 1.0])
}

#[test]
fn test_cr_values() {
    let s = simple_system("test");
    // CR = λ₂/λₙ = 0.5/5.0 = 0.1
    assert!((s.cr() - 0.1).abs() < 1e-10);
}

#[test]
fn test_spectral_gap() {
    let s = simple_system("test");
    assert!((s.spectral_gap() - 0.5).abs() < 1e-10);
}

#[test]
fn test_mode_energy() {
    let s = SpectralSystem::new("test", vec![0.0, 1.0], vec![3.0, 4.0]);
    assert!((s.mode_energy() - 25.0).abs() < 1e-10);
}

#[test]
fn test_cosine_similarity_same() {
    assert!((cosine_similarity(&[1.0, 2.0], &[1.0, 2.0]) - 1.0).abs() < 1e-10);
}

#[test]
fn test_cosine_similarity_orthogonal() {
    assert!(cosine_similarity(&[1.0, 0.0], &[0.0, 1.0]).abs() < 1e-10);
}

#[test]
fn test_cosine_similarity_empty() {
    assert_eq!(cosine_similarity(&[], &[]), 0.0);
}

#[test]
fn test_couple_reinforcing() {
    let a = simple_system("A");
    let b = simple_system("B");
    let result = couple(&a, &b, 1.0);
    assert_eq!(result.mode, CouplingMode::Reinforcing);
    assert!(result.alignment > 0.9);
}

#[test]
fn test_couple_competitive() {
    let a = SpectralSystem::new("A", vec![0.0, 0.5, 2.0], vec![1.0, 1.0, 1.0]);
    let b = SpectralSystem::new("B", vec![0.0, 0.5, 2.0], vec![-1.0, -1.0, -1.0]);
    let result = couple(&a, &b, 1.0);
    assert!(result.alignment < -0.9);
}

#[test]
fn test_emergence_curve_length() {
    let a = simple_system("A");
    let b = simple_system("B");
    let curve = emergence_curve(&a, &b, 10);
    assert_eq!(curve.len(), 11); // 0..=10
}

#[test]
fn test_critical_coupling() {
    let w = SpectralSystem::new("W", vec![0.0, 0.5, 2.0, 5.0], vec![1.0, 0.0]);
    let crit = critical_coupling(&w, &w, 50);
    assert!(crit.is_some());
    let val = crit.unwrap();
    assert!(val >= 0.0 && val <= 1.0);
}

#[test]
fn test_graph_to_system() {
    let adj = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]];
    let sys = graph_to_system("path", &adj);
    assert_eq!(sys.eigenvalues.len(), 3);
    assert!(sys.eigenvalues[0].abs() < 1e-8, "smallest eigenvalue should be ~0");
}

#[test]
fn test_chain_coupling() {
    let systems = vec![
        simple_system("A"),
        simple_system("B"),
        SpectralSystem::new("C", vec![0.0, 0.8, 3.0, 6.0], vec![1.0, -1.0, 1.0, -1.0]),
    ];
    let result = chain_coupling(&systems, 0.5);
    assert!(!result.coupled.eigenvalues.is_empty());
    assert!(result.synergy > 0.0);
}
