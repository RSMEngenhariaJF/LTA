# Circuitos com linhas convencionais — Roteiro (Pozar)

Fonte: D. M. Pozar, *Microwave Engineering*, 4ª ed., Wiley, 2012. Complementa o capítulo 0 (linhas convencionais) e as aplicações do capítulo 4 do Martín.

| Seção | Tema | Situação no LTA (v1.9) |
|---|---|---|
| 2.4 | Carta de Smith | Na aba Casamento (Γin(f), carga e solução) e nos Ressonadores |
| 3.6, 3.8 | Ondas de superfície, modos superiores, fórmulas da microstrip | Limites TE1/TM1 e modo transversal com aviso; eqs. 3.195–3.199 para comparação (aba Substrato) |
| 4.6 | Descontinuidades | Ponta aberta (Kirschning-Jansen-Koster) e gap série (Garg-Bahl); junção T, curva e degrau pendentes |
| 5.1 | Rede L | Ex. 5.1, duas topologias e componentes em f0 |
| 5.2–5.3 | Stub simples (derivação e série) e stub duplo | Ex. 5.2–5.4, com layout e desconto da ponta aberta |
| 5.4–5.7 | Transformadores λ/4, binomial e Chebyshev | Ex. 5.6–5.7 e Tabela 5.2; N qualquer |
| 5.8 | Tapers exponencial, triangular e Klopfenstein | Ex. 5.8; resposta exata por cascata; layout contínuo |
| 5.9 | Bode-Fano | RC e RL série |
| 6.2 | Ressonadores de linha, Q0 | Ex. 6.2 (o Q0 do livro herda um β inconsistente; ver manual) |
| 6.6 | Ressonador acoplado por gap | Ex. 6.6; gap físico e layout |
| 7.2–7.3 | Junção T e Wilkinson (igual e desigual) | Ex. 7.1 e eq. 7.37; transformadores de saída opcionais; layouts |
| 7.5, 7.8 | Branch-line e anel híbrido | Microstrip convencional e slow-wave (LTA); layouts |
| 7.6 | Acopladores de linhas acopladas, 1 e N seções | Ex. 7.7–7.8; versão LTA com modo ímpar compensado |
| 7.7 | Lange | Eqs. 7.97–7.100; sem layout (fios de solda) |
| 8.1–8.2 | Estruturas periódicas e parâmetros imagem | Cobertos pelo cap. 2 do Martín (Bloch); método imagem pendente |
| 8.3–8.5 | Perda de inserção, Richards e Kuroda | Ex. 8.5 (N ≤ 3) |
| 8.6 | Degraus de impedância | Ex. 8.6, com protótipo concentrado para comparação |
| 8.7 | Passa-faixa de linhas acopladas | Ex. 8.7; versão LTA suprime a espúria em 2f0 |
| 8.8 | Stubs λ/4 e ressonadores acoplados por capacitores | Ex. 8.8–8.10 |
| 9–14 | Ferrites, ruído, ativos, amplificadores, osciladores, sistemas | Fora do escopo |

## Combinações com linhas artificiais
- **Acoplador de linhas acopladas LTA:** capacitores periódicos entre as linhas carregam só o modo ímpar e igualam as velocidades de fase; a diretividade da microstrip sobe de ~12 para ~70 dB (10 dB, RO4003C).
- **Filtro de linhas acopladas LTA:** a mesma compensação remove a banda espúria em 2f0 (Martín, seção 4.2.3.2).
- **Híbridas e divisores slow-wave:** ramos λ/4 trocados por linhas com carga capacitiva periódica (stubs ou patches), com fase e impedância de Bloch exatas em f0; ramos 40–55% mais curtos.
- As implementações CRLH dual-band e LH de banda larga continuam disponíveis e agora são comparadas com a versão em microstrip convencional.

## Notas
- Ex. 6.2: β = 151,0 rad/m no livro; com os dados do exemplo, β = 140,5 rad/m e Q0 ≈ 729.
- Ex. 7.8: a teoria de acoplamento fraco dá 20 dB em f0; a cascata exata, 19,97 dB.
- Garg-Bahl (microstrip acoplada) tem precisão de ~3% em 0,2 ≤ W/h ≤ 2 e 0,05 ≤ S/h ≤ 2.
