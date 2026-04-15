# Percolação em Rede 2D + SIR + Fractal

## Arquivos

**`percolacao.py`** — toda a lógica computacional. Sem dependências externas.

**`index.html`** — interface visual. Carrega `percolacao.py` via Pyodide (Python no navegador). Não possui lógica. 

**`Percolação_e_dimensão_fractal.pdf`** - relatório com toda descrição e resultados numéricos.

## Como rodar

```bash
python3 -m http.server 8000
# abrir: http://localhost:8000/index.html
```

O servidor é necessário para o `fetch('percolacao.py')` funcionar.

## Funções

**`make_grid(n, p)`** — grade n×n onde cada célula é ocupada com probabilidade p.

**`label_clusters(grid, n)`** — Union-Find com path compression e union by rank. Retorna o rótulo canônico de cada célula e o tamanho de cada cluster. Complexidade O(N²).

**`find_perc_cluster(labels, sizes, n)`** — verifica se algum cluster toca simultaneamente o topo e a base da grade.

**`box_counting(labels, perc_label, n)`** — cobre o cluster percolante com caixas de lado 2^b para b = 0, 1, … ,log₂(n). Retorna pares (eps, N).

**`fractal_dimension(box_data)`** — regressão linear em ln(N) × ln(1/eps). A inclinação é Df.

**`init_sir(grid, n)`** — converte a grade de percolação em grade SIR (0=imune, 1=suscetível).

**`step_sir(sir_grid, n)`** — um passo do motor SIR: cada infectado (2) contamina os quatro vizinhos suscetíveis e passa a recuperado (3) no mesmo passo.

**`count_sir_states(sir_grid)`** — conta S, I, R.

**`perc_threshold_day(history_s, reachable_s0)`** — retorna o dia em que 59,27% dos suscetíveis alcançáveis foram removidos (limiar p_c = 0.592746).

## Validar

```bash
python3 percolacao.py
```