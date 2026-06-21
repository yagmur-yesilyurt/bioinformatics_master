"""
DBSCAN CLUSTERING ANALYSIS: Conformational States
Identifies distinct conformational clusters in MD trajectory using PCA and DBSCAN

OUTPUT FILE:
------------
conformational_clustering_dbscan.png - PCA scatter plot with identified clusters
"""

import MDAnalysis as mda
from MDAnalysis.analysis.pca import PCA
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import warnings


# ============================================
# CONFIGURATION
# ============================================
class Config:
    """Analysis parameters"""
    PSF_FILE = "step5_input.psf"
    N_FILES = 10
    SELECTION = "name CA"  # Analyze C-alpha atoms only
    
    # PCA parameters
    N_COMPONENTS = 2
    ALIGN = True
    
    # DBSCAN parameters
    # eps: Maximum distance for points to be neighbors
    # Higher values = larger clusters, lower values = more noise
    EPS = 10
    
    # min_samples: Minimum points to form a cluster
    # Higher values = denser clusters required
    MIN_SAMPLES = 50


# ============================================
# UTILITY FUNCTIONS
# ============================================
def suppress_warnings():
    """Suppress unnecessary warnings"""
    warnings.filterwarnings('ignore')


def generate_trajectory_list(n_files):
    """Generate list of DCD trajectory files"""
    return [f"step7_{i}.dcd" for i in range(1, n_files + 1)]


def print_section(title):
    """Print formatted section header"""
    print(f"\n{'='*60}")
    print(title)
    print('='*60)


# ============================================
# DATA LOADING AND PCA
# ============================================
def load_and_transform(config):
    """Load trajectory and perform PCA transformation"""
    print(">>> Loading data and preparing PCA...")
    
    # Generate trajectory file list
    dcd_files = generate_trajectory_list(config.N_FILES)
    
    # Load universe
    u = mda.Universe(config.PSF_FILE, dcd_files)
    print(f"✓ Loaded {len(u.trajectory):,} frames from {config.N_FILES} files")
    
    # Perform PCA
    print(f">>> Running PCA on {config.SELECTION} atoms...")
    pc = PCA(
        u,
        select=config.SELECTION,
        align=config.ALIGN,
        mean=None,
        n_components=config.N_COMPONENTS
    ).run()
    
    # Transform coordinates to PC space
    pca_space = pc.transform(
        u.select_atoms(config.SELECTION),
        n_components=config.N_COMPONENTS
    )
    
    print(f"✓ PCA complete")
    print(f"  PC1 variance: {pc.cumulated_variance[0]*100:.1f}%")
    print(f"  PC2 variance: {(pc.cumulated_variance[1]-pc.cumulated_variance[0])*100:.1f}%")
    
    return pca_space, pc


# ============================================
# CLUSTERING ANALYSIS
# ============================================
def perform_clustering(pca_space, config):
    """Apply DBSCAN clustering to PCA-transformed data"""
    print(f"\n>>> Running DBSCAN clustering...")
    print(f"    Parameters: eps={config.EPS}, min_samples={config.MIN_SAMPLES}")
    
    db = DBSCAN(eps=config.EPS, min_samples=config.MIN_SAMPLES).fit(pca_space)
    labels = db.labels_
    
    # Calculate statistics
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    n_total = len(labels)
    
    print(f"\n✓ Clustering complete:")
    print(f"  Identified clusters: {n_clusters}")
    print(f"  Noise points: {n_noise:,} ({n_noise/n_total*100:.1f}%)")
    print(f"  Clustered points: {n_total-n_noise:,} ({(n_total-n_noise)/n_total*100:.1f}%)")
    
    return labels, n_clusters, n_noise


# ============================================
# VISUALIZATION
# ============================================
def plot_clusters(pca_space, labels, pc, n_clusters, n_noise, config, output_file):
    """Create publication-quality cluster visualization"""
    print("\n>>> Creating visualization...")
    
    plt.figure(figsize=(10, 8))
    
    # Generate colors for each cluster
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) 
              for each in np.linspace(0, 1, len(unique_labels))]
    
    # Plot each cluster
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Noise points: black, semi-transparent, small
            col = [0, 0, 0, 1]
            label_name = "Transition/Noise"
            marker_size = 3
            alpha_val = 0.3
        else:
            # Cluster points: colored, opaque, larger
            label_name = f"Cluster {k}"
            marker_size = 6
            alpha_val = 0.8
        
        class_member_mask = (labels == k)
        xy = pca_space[class_member_mask]
        
        plt.plot(xy[:, 0], xy[:, 1], 'o',
                markerfacecolor=tuple(col),
                markeredgecolor='none',
                markersize=marker_size,
                alpha=alpha_val,
                label=label_name)
    
    # Labels and formatting
    plt.title(
        f'DBSCAN Clustering (eps={config.EPS})\n'
        f'Clusters: {n_clusters}, Noise Points: {n_noise:,}',
        fontsize=14, fontweight='bold'
    )
    plt.xlabel(f'PC1 ({pc.cumulated_variance[0]*100:.1f}%)', fontsize=12)
    plt.ylabel(
        f'PC2 ({(pc.cumulated_variance[1]-pc.cumulated_variance[0])*100:.1f}%)',
        fontsize=12
    )
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='best', fontsize=9)
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Figure saved: {output_file}")


# ============================================
# MAIN EXECUTION
# ============================================
def main():
    """Main analysis pipeline"""
    print_section("DBSCAN CLUSTERING ANALYSIS")
    
    # Initialize configuration
    config = Config()
    suppress_warnings()
    
    # Load data and perform PCA
    pca_space, pc = load_and_transform(config)
    
    # Perform clustering
    labels, n_clusters, n_noise = perform_clustering(pca_space, config)
    
    # Create visualization
    output_file = "08_conformational_clustering_dbscan.png"
    plot_clusters(pca_space, labels, pc, n_clusters, n_noise, config, output_file)
    
    print_section("ANALYSIS COMPLETE")
    print(f"\nResults saved to: {output_file}")
    print("\nClustering Summary:")
    print(f"  • Total conformations: {len(labels):,}")
    print(f"  • Distinct states: {n_clusters}")
    print(f"  • Transition frames: {n_noise:,}")


if __name__ == "__main__":
    main()
