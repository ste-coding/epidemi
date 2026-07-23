# Sarampo: SEIR de campo médio vs. autômato celular espacial, ligados pela dimensão de similaridade

**Entrega conjunta para duas disciplinas na UFRPE (BSI):**
- Modelagem Matemático-Computacional Aplicada à Epidemiologia
- Autômatos Celulares

---

## O que é isto

A mesma doença (**sarampo**) é modelada de duas formas, com **parâmetros
epidemiológicos idênticos**, e as duas são conectadas por um conceito da teoria
de autômatos celulares: a **dimensão de similaridade (fractal)**.

| | Modelo | Disciplina | Referência |
|---|---|---|---|
| **Contínuo** | **SEIR** de campo médio (EDOs, Runge–Kutta) | Epidemiologia | Keeling & Rohani, Cap. 1–2 |
| **Discreto/espacial** | **Autômato celular SEIR** 2D | Autômatos Celulares | Fuentes & Kuperman (1999); Schiff, Cap. 1 |
| **Ponte** | **Dimensão de similaridade** do padrão epidêmico (contagem de caixas, validada na curva de Koch) | ambas | Schiff, Cap. 1 (Eq. `D_s = log N / log(1/r)`) |

### A história científica
Com a *mesma* transmissibilidade por contato, restringir o contato a uma
vizinhança local converte o surto abrupto de campo médio (pico de prevalência
**25,4%** no dia 20) em uma onda viajante lenta (pico de **1,5%** no dia 494).
Uma varredura de parâmetros revela um **limiar de invasão do tipo percolação**
(τ\* ≈ 0,07) que o limiar suave `R0 = 1` da EDO nunca mostra. A **dimensão de
similaridade** quantifica a geometria: a frente ativa é quase 1D (`D ≈ 1,1`),
enquanto o aglomerado de já infectados cresce de um fractal (`D ≈ 1,3`, perto do
limiar) rumo ao valor que preenche o plano à medida que a transmissibilidade
aumenta, passando pelo valor de percolação 2D `≈ 1,90`.

## Passos para reprodução

```bash
pip install -r requirements.txt
cd src
python run_experiments.py          # gera figures/ e results.json
```

Cada módulo também roda isolado, para uma checagem rápida:

```bash
python src/seir_ode.py
python src/seir_ca.py
python src/similarity_dimension.py   # imprime Koch D_s = log4/log3
```

## Parâmetros principais (sarampo)

| símbolo | significado | valor |
|---|---|---|
| `R0` | número básico de reprodução | 15 |
| `1/gamma` | período infeccioso | 5 dias |
| `1/sigma` | período latente | 8 dias |
| `beta = R0·gamma` | taxa de transmissão | 3,0 /dia |
| `tau` | prob. de transmissão por contato (AC) | varrida de 0,03 a 0,70 |

## Referências
1. M. J. Keeling & P. Rohani. *Modeling Infectious Diseases in Humans and Animals.* Princeton Univ. Press, 2008 (Cap. 1–2).
2. J. L. Schiff. *Cellular Automata: A Discrete View of the World.* Wiley, 2008 (Cap. 1).
3. M. A. Fuentes & M. N. Kuperman. Cellular automata and epidemiological models with spatial dependence. *Physica A* **267**, 471–486 (1999).
4. W. O. Kermack & A. G. McKendrick. A contribution to the mathematical theory of epidemics. *Proc. R. Soc. A* **115**, 700–721 (1927).
