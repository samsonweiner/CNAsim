import numpy as np

def choose_noise_model(choice, tree, chrom_names, err_rate, dub_rate):
    if choice == 1:
        add_noise_simple(tree, chrom_names, err_rate)
    elif choice == 2:
        add_noise_normal(tree, chrom_names, err_rate)
    elif choice == 3:
        add_noise_boundary(tree, chrom_names, err_rate)
    elif choice == 4:
        add_noise_sequence(tree, chrom_names, err_rate, dub_rate)
    elif choice == 5:
        add_noise_boundary(tree, chrom_names, err_rate)
        add_noise_normal(tree, chrom_names, dub_rate)

def add_noise_simple(tree, chrom_names, err_rate):
    for leaf in tree.iter_leaves():
        for chrom in chrom_names:
            for allele in [0, 1]:
                for b in range(len(leaf.profile[chrom][allele])):
                    if np.random.uniform() < err_rate:
                        if leaf.profile[chrom][allele][b] != 0:
                            if np.random.binomial(1, 0.5) == 0:
                                leaf.profile[chrom][allele][b] += 1
                            else:
                                leaf.profile[chrom][allele][b] -= 1

def add_noise_normal(leaf, chrom_names, err_rate):
    for chrom in chrom_names:
        for allele in [0, 1]:
            for b in range(len(leaf.profile[chrom][allele])):
                leaf.profile[chrom][allele][b] = round(np.random.normal(leaf.profile[chrom][allele][b], err_rate*leaf.profile[chrom][allele][b]))

def add_noise_boundary(leaf, chrom_names, err_rate):
    for chrom in chrom_names:
        for allele in [0, 1]:
            profile = leaf.profile[chrom][allele]
            i, seg_length = 1, 1
            prev_num = profile[0]
            while i < len(profile):
                if profile[i] != prev_num:
                    new_length = max(round(np.random.normal(seg_length, seg_length*err_rate)), 1)
                    if new_length <= seg_length:
                        for j in range(i - seg_length + new_length, i):
                            profile[j] = profile[i]
                        seg_length = 1
                        prev_num = profile[i]
                        i += 1
                    else:
                        for j in range(i, min(i + new_length - seg_length, len(profile))):
                            profile[j] = prev_num
                        seg_length = 1
                        i += new_length - seg_length
                else:
                    seg_length += 1
                    i += 1

def add_noise_sequence(tree, chrom_names, cov_rate, dub_rate):
    #select doublets
    doublets = {}
    cell_names = list(tree.iter_leaves())
    for i in range(len(cell_names)):
        if np.random.uniform() < dub_rate:
            pooled_cell = cell_names[np.random.choice([j for j in range(len(cell_names)) if j != i])]
            doublets[cell_names[i]] = pooled_cell
    for chrom in chrom_names:
        for allele in [0, 1]:
            for doublet, pooled_cell in doublets.items():
                profile = doublet.profile[chrom][allele]
                alt_profile = pooled_cell.profile[chrom][allele]
                for i in range(len(profile)):
                    if np.random.binomial(1, 0.5) == 0:
                        profile[i] = alt_profile[i]
                    
    #coverage non-uniformity
    for leaf in tree.iter_leaves():
        for chrom in chrom_names:
            for allele in [0, 1]:
                profile = leaf.profile[chrom][allele]
                cov_chart = gen_coverage(len(profile), np.random.randint(4, 7))
                for i in range(len(profile)):
                    copy_num = profile[i]
                    coverage = 2*cov_chart[i]
                    if coverage < 1:
                        scale = 1 - (coverage*cov_rate)
                        profile[i] = round(copy_num*scale)
                    elif coverage > 1:
                        scale = 1 + (coverage*cov_rate)
                        profile[i] = round(copy_num*scale)

# Jitter and boundary model combined
def add_noise_mixed(tree, chrom_names, boundary_err_rate, jitter_err_rate):
    for leaf in tree.iter_leaves():
        if leaf.cell_type != 'normal':
            add_noise_boundary(leaf, chrom_names, boundary_err_rate)
            add_noise_normal(leaf, chrom_names, jitter_err_rate)


def bezier_coef(points):
    n = len(points) - 1

    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2

    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]

    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2

    return A, B

def single_cubic_bezier(a, b, c, d):
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d

def multiple_cubic_bezier(points):
    A, B = bezier_coef(points)
    return [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]

def gen_random_interval(num_bins, scale):
    x = [i for i in range(0, num_bins, round(num_bins/scale))]
    if x[-1] != num_bins-1:
        if num_bins - x[-1] >= round(num_bins/scale)/2:
            x.append(num_bins-1)
        else:
            x[-1] = num_bins-1
    y = [np.random.uniform(0, 1) for i in range(len(x))]
    points = [np.array(p) for p in zip(x, y)]
    return points
    
def gen_coverage(num_bins, scale):
    points = gen_random_interval(num_bins, scale)
    A, B = bezier_coef(points)
    curves = [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]

    new_points = []
    for i in range(len(points) - 1):
        f = curves[i]
        gaps = points[i+1][0] - points[i][0] + 1
        coords = [f(t) for t in np.linspace(0, 1, int(gaps))]
        if i == 0:
            new_points.append((round(coords[0][0]), coords[0][1]))
        new_points += [(round(i[0]), i[1]) for i in coords[1:]]
    
    new_points.sort(key=lambda x: x[0])
    new_points = [max(min(i[1], 1), 0) for i in new_points]
    return new_points
