import mdtraj as md
import numpy as np
from sklearn.cluster import DBSCAN
from collections import Counter

def get_most_probable_conformation(
    trajectory_file,
    topology_file,
    output_path,
    selection='protein and backbone',
    eps=1.0,
    min_samples=10,
    subsample_rate=1,
    # --- YENİ PARAMETRELER ---
    ligand_resname='UNL',           # Ligand residue ismi
    ploop_resids=[348, 349, 353, 355],  # PRO348, ASP349, THR353, PHE355
    ligand_cutoff_angstrom=10.0     # Angstrom cinsinden filtreleme mesafesi
):
    print(f"Loading trajectory: {trajectory_file}")
    t = md.load(trajectory_file, top=topology_file)
    print(f"Total frames loaded: {len(t)}")

    # ─────────────────────────────────────────────────────────────
    # ADIM 1: Ligand P-loop mesafesine göre frame filtreleme
    # ─────────────────────────────────────────────────────────────
    ligand_cutoff_nm = ligand_cutoff_angstrom / 10.0  # mdtraj nm kullanır

    lig_indices = t.top.select(f'resname {ligand_resname}')
    ploop_sel = ' or '.join([f'resSeq {r}' for r in ploop_resids])
    ploop_indices = t.top.select(ploop_sel)

    if len(lig_indices) == 0 or len(ploop_indices) == 0:
        print("WARNING: Ligand veya P-loop atomları bulunamadı — filtreleme atlanıyor.")
        filtered_t = t
        frame_map = np.arange(len(t))
    else:
        print(f"Ligand atomu sayısı: {len(lig_indices)}")
        print(f"P-loop atom sayısı: {len(ploop_indices)} (resids: {ploop_resids})")
        print(f"Filtreleme: mexiletine < {ligand_cutoff_angstrom:.1f} Å olan frame'ler...")

        # Tüm ligand-ploop çiftleri
        pairs = np.array([[l, p] for l in lig_indices for p in ploop_indices])

        # Bellek dostu: 500'er frame'lik parçalarda hesapla
        CHUNK = 500
        min_dists = np.zeros(len(t))
        for start in range(0, len(t), CHUNK):
            end = min(start + CHUNK, len(t))
            chunk_dists = md.compute_distances(t[start:end], pairs)
            min_dists[start:end] = chunk_dists.min(axis=1)

        mask = min_dists < ligand_cutoff_nm
        n_filtered = mask.sum()
        print(f"Filtreyi geçen frame sayısı: {n_filtered} / {len(t)} ({100*n_filtered/len(t):.1f}%)")

        if n_filtered == 0:
            print(f"UYARI: Hiçbir frame filtreyi geçemedi!")
            print(f"Tüm frame'lerdeki min mesafe: {min_dists.min()*10:.2f} Å")
            print("Filtreleme atlanıyor — cutoff'u artırmayı dene.")
            filtered_t = t
            frame_map = np.arange(len(t))
        else:
            frame_map = np.where(mask)[0]
            filtered_t = t[frame_map]

    # ─────────────────────────────────────────────────────────────
    # ADIM 2: Backbone RMSD üzerinden DBSCAN clustering
    # ─────────────────────────────────────────────────────────────
    selection_indices = filtered_t.top.select(selection)
    t_selection = filtered_t.atom_slice(selection_indices)
    sub_t = t_selection[::subsample_rate]
    sub_frame_map = frame_map[::subsample_rate]
    n_frames = len(sub_t)
    print(f"\nClustering için frame sayısı: {n_frames}")

    # RMSD matrisi
    rmsd_matrix = np.zeros((n_frames, n_frames))
    for i in range(n_frames):
        if i % 500 == 0:
            print(f"RMSD Matrix: {i}/{n_frames}...")
        rmsd_matrix[i] = md.rmsd(sub_t, sub_t, frame=i)

    print("DBSCAN clustering yapılıyor...")
    db = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed', n_jobs=-1).fit(rmsd_matrix)
    labels = db.labels_

    cluster_counts = Counter(labels)
    clean_counts = {k: v for k, v in cluster_counts.items() if k != -1}

    print(f"\nCluster sayısı: {len(clean_counts)}")
    print(f"Noise frame sayısı: {cluster_counts.get(-1, 0)}")
    for label, count in sorted(clean_counts.items(), key=lambda x: -x[1]):
        print(f"  Cluster {label}: {count} frame ({100*count/n_frames:.1f}%)")

    if not clean_counts:
        print("HATA: Hiç cluster bulunamadı. eps değerini artır.")
        return None

    # En büyük cluster'dan medoid seç (diğer frame'lere ortalama RMSD'si en düşük frame)
    largest_cluster_label = max(clean_counts, key=clean_counts.get)
    indices_in_cluster = np.where(labels == largest_cluster_label)[0]

    cluster_rmsd = rmsd_matrix[np.ix_(indices_in_cluster, indices_in_cluster)]
    avg_rmsd = cluster_rmsd.mean(axis=1)
    medoid_idx = np.argmin(avg_rmsd)
    rep_idx_in_sub = indices_in_cluster[medoid_idx]

    print(f'Medoid ortalama RMSD: {avg_rmsd[medoid_idx]*10:.3f} Å')

    # Orijinal trajectory index'e geri dön
    rep_idx_original = sub_frame_map[rep_idx_in_sub]

    # Tam atomlarla kaydet
    output_file = f'{output_path}/rep_structure_full_10k.pdb'
    t[rep_idx_original].save(output_file)

    print(f'\nTemsilci frame kaydedildi: {output_file}')
    print(f'Orijinal trajectory index: {rep_idx_original}')
    print(f'Cluster boyutu: {clean_counts[largest_cluster_label]} / {n_frames}')
    return rep_idx_original
