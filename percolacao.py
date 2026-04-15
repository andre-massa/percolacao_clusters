"""
percolacao.py
=============
Funções computacionais do modelo de percolação + simulação SIR.

Não há nenhuma dependência de interface gráfica aqui.

Constante central:
    PC_SITE_SQUARE = 0.592746  — limiar crítico de percolação de sítios
    na rede quadrada 2D, determinado por simulações de Monte Carlo de
    alta precisão (Newman & Ziff, 2000: pc = 0.59274602 ± 0.00000004).
"""

import random
import math

PC_SITE_SQUARE = 0.592746


def make_grid(n, p):
    """
    Gera uma grade nxn onde cada célula é ocupada (1) com probabilidade p.
    Retorna uma lista plana de inteiros (0 ou 1) com n*n elementos.
    Índice: i = row * n + col
    """
    return [1 if random.random() < p else 0 for _ in range(n * n)]


def label_clusters(grid, n):
    """
    Encontra todos os clusters de células ocupadas usando Union-Find.

    Retorna:
        labels : lista com o rótulo (raiz canônica) de cada célula,
                 -1 para células vazias.
        sizes  : dicionário {raiz: tamanho} de cada cluster.
    """
    parent = list(range(n * n))   # cada célula começa sendo sua própria raiz
    rank   = [0] * (n * n)

    def find(x):
        # path compression
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a == b:
            return
        # union by rank
        if rank[a] < rank[b]:
            a, b = b, a
        parent[b] = a
        if rank[a] == rank[b]:
            rank[a] += 1

    # une vizinhos ocupados (4-conectividade: cima e esquerda)
    for i in range(n * n):
        if not grid[i]:
            continue
        r, c = divmod(i, n)
        if r > 0 and grid[i - n]:
            union(i, i - n)
        if c > 0 and grid[i - 1]:
            union(i, i - 1)

    # monta labels e conta tamanhos
    labels = [-1] * (n * n)
    sizes  = {}
    for i in range(n * n):
        if grid[i]:
            root = find(i)
            labels[i] = root
            sizes[root] = sizes.get(root, 0) + 1

    return labels, sizes


def find_perc_cluster(labels, sizes, n):
    """
    Verifica se existe um cluster que toca a linha 0 (topo)
    e a linha n-1 (base) ao mesmo tempo.

    Retorna o rótulo do maior cluster percolante, ou -1 se não existir.
    """
    top_roots = set()
    for c in range(n):
        lbl = labels[c]          # linha 0
        if lbl >= 0:
            top_roots.add(lbl)

    best, best_size = -1, 0
    for c in range(n):
        lbl = labels[(n - 1) * n + c]   # linha n-1
        if lbl >= 0 and lbl in top_roots:
            s = sizes.get(lbl, 0)
            if s > best_size:
                best_size = s
                best = lbl

    return best



def box_counting(labels, perc_label, n):
    """
    Estima a dimensão fractal do cluster de percolação pelo método de box-counting.

    Para cada escala b = 1, 2, ..., log2(n):
        - cobre o cluster com caixas de lado 2^b
        - conta quantas caixas contêm ao menos uma célula do cluster

    Retorna lista de dicionários: [{'eps': float, 'N': int}, ...]
        eps = tamanho relativo da caixa (2^b / n)
        N   = número de caixas ocupadas
    """
    if perc_label < 0:
        return []

    cells = [(i // n, i % n) for i, lbl in enumerate(labels) if lbl == perc_label]
    if not cells:
        return []

    results = []
    max_pow = int(math.log2(n))

    for b in range(0, max_pow + 1):
        bs = 2 ** b
        boxes = set()
        for r, c in cells:
            boxes.add((r // bs, c // bs))
        results.append({'eps': bs / n, 'N': len(boxes)})

    return results


def lin_reg(xs, ys):
    """
    Regressão linear simples por mínimos quadrados.
    Retorna (slope, intercept).
    """
    n   = len(xs)
    sx  = sum(xs)
    sy  = sum(ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    sx2 = sum(x * x for x in xs)

    slope     = (n * sxy - sx * sy) / (n * sx2 - sx * sx)
    intercept = (sy - slope * sx) / n
    return slope, intercept


def fractal_dimension(box_data):
    """
    Calcula Df a partir dos dados de box-counting.
    Df = inclinação de ln(N) x ln(1/eps).
    """
    if len(box_data) < 4:
        data = box_data
    else: 
        data = box_data[1:-1]  
        
    if len(data) < 2:
        return float('nan')  # não dá para estimar Df com menos de 2 pontos
    xs = [math.log(1 / d['eps']) for d in box_data]
    ys = [math.log(d['N'])       for d in box_data]
    slope, _ = lin_reg(xs, ys)
    return slope



def init_sir(grid, n):
    """
    Inicializa a simulação SIR a partir da grade de percolação.

    Estados:
        0 = imune - (célula vazia) 
        1 = suscetível
        2 = infectado
        3 = recuperado

    Retorna a grade SIR como lista de inteiros.
    """
    return [1 if cell else 0 for cell in grid]


def step_sir(sir_grid, n):
    """
    Executa um passo da simulação SIR determinística:
        - Cada infectado (2) infecta todos os vizinhos suscetíveis (1) em cruz (4-conectividade)
        - Depois passa imediatamente para recuperado (3)

    Retorna:
        new_grid      : nova grade SIR após o passo
        newly_infected: número de novas infecções neste passo
        has_active    : True se ainda há infectados ativos
    """
    new_grid = list(sir_grid)
    newly_infected = 0
    has_active = False

    for i in range(n * n):
        if sir_grid[i] != 2:
            continue
        has_active = True
        r, c = divmod(i, n)

        # infecta vizinhos suscetíveis (4-conectividade)
        neighbors = []
        if r > 0:     neighbors.append(i - n)
        if r < n - 1: neighbors.append(i + n)
        if c > 0:     neighbors.append(i - 1)
        if c < n - 1: neighbors.append(i + 1)

        for nb in neighbors:
            if sir_grid[nb] == 1 and new_grid[nb] == 1:
                new_grid[nb] = 2
                newly_infected += 1

        # infectado passa a recuperado no mesmo passo
        new_grid[i] = 3

    return new_grid, newly_infected, has_active


def perc_threshold_day(history_s, reachable_s0):
    """
    Encontra o dia em que a simulação SIR cruza o limiar crítico pc.

    O critério é: quando a fração de suscetíveis ALCANÇÁVEIS já removidos
    (infectados + recuperados) ultrapassa PC_SITE_SQUARE = 0.592746.

    Por que usar reachable_s0 e não o total da grade?
    Suscetíveis em clusters isolados nunca são atingidos pela infecção.
    Usá-los como base inflaria o denominador e o threshold nunca seria
    cruzado. A base correta é o tamanho do cluster onde a infecção corre.

    Parâmetros:
        history_s    : lista com S(t) a cada passo da simulação
        reachable_s0 : tamanho do cluster alcançável (candidatos ao paciente zero)

    Retorna:
        Dia (índice inteiro) em que pc é cruzado, ou -1 se nunca ocorrer.
    """
    if not history_s or reachable_s0 <= 0:
        return -1

    threshold = reachable_s0 * PC_SITE_SQUARE

    # Desconto de suscetíveis isolados (não alcançáveis)
    isolated = history_s[0] - reachable_s0

    for day, s in enumerate(history_s[1:], start=1):
        removed = reachable_s0 - (s - isolated)
        if removed >= threshold:
            return day

    return -1


def count_sir_states(sir_grid):
    """
    Conta as células em cada estado SIR.
    Retorna (suscetíveis, infectados, recuperados).
    """
    s = sir_grid.count(1)
    i = sir_grid.count(2)
    r = sir_grid.count(3)
    return s, i, r


if __name__ == "__main__":
    import time

    print(f"pc (rede quadrada 2D, sítios) = {PC_SITE_SQUARE}")
    print()

    for N, p in [(100, 0.59), (200, 0.593), (300, 0.592746)]:
        t0 = time.perf_counter()
        grid = make_grid(N, p)
        labels, sizes = label_clusters(grid, N)
        perc = find_perc_cluster(labels, sizes, N)
        elapsed = (time.perf_counter() - t0) * 1000

        percolou = "SIM" if perc >= 0 else "NÃO"
        n_clusters = len(sizes)
        maior = max(sizes.values()) if sizes else 0
        print(f"N={N:3d}  p={p:.6f}  percolou={percolou}  "
              f"clusters={n_clusters:5d}  maior={maior:6d}  {elapsed:.1f}ms")

        if perc >= 0:
            bc = box_counting(labels, perc, N)
            df = fractal_dimension(bc)
            print(f"         Df = {df:.4f}  (esperado ≈ 1.89 em pc)")

        # simula SIR e encontra dia do pc
        sir = init_sir(grid, N)
        candidates = [i for i, lbl in enumerate(labels) if lbl == perc] if perc >= 0 else []
        if candidates:
            idx = candidates[len(candidates) // 2]
            sir[idx] = 2
            hist_s = [sir.count(1)]
            hist_i = [1]
            for _ in range(500):
                sir, new_inf, active = step_sir(sir, N)
                s, i, r = count_sir_states(sir)
                hist_s.append(s)
                hist_i.append(i)
                if not active or i == 0:
                    break
            pc_day = perc_threshold_day(hist_s, len(candidates))
            pico   = hist_i.index(max(hist_i))
            print(f"         SIR: {len(hist_s)-1} dias  "
                  f"pico=dia {pico}  pc_day=dia {pc_day}")
        print()