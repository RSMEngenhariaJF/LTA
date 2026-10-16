# Linhas de Transmissão Metamateriais — Roteiro (capítulo 3)

Fonte: F. Martín, *Artificial Transmission Lines for RF and Microwave Applications*, Wiley, 2015, cap. 3 e Apêndices G e H.

> Documento de trabalho. Equações conferidas contra o PDF; trechos com texto ambíguo estão anotados.

| Seção | Tema | Situação no LTA (v1.4) |
|---|---|---|
| 3.2 | Meios efetivos LH (fio + SRR) | Não implementado (sem aplicação direta em projeto de linhas) |
| 3.3 | Ressonadores pequenos (SRR, CSRR, OSRR, OCSRR…) | CSRR: estimativa quase-estática de fo (dualidade) e desenho das fendas |
| 3.4.1–3.4.2 | Linha dual e CRLH | Aba CRLH: eqs. 3.51–3.66, síntese balanceada por banda e por fases, meio efetivo, leaky-wave |
| 3.4.3 | D-CRLH e ordens superiores | D-CRLH implementada; ordens 3 e 4 pendentes |
| 3.5.1 | CL-loaded (interdigital + stub com via) | Aba CRLH microstrip CL-loaded: geometria inicial, layout com vias, extração |
| 3.5.2.1 | SRR em CPW | Extração (SRR e SRR + fita) e transformação 3.71–3.74; projeto geométrico pendente (falta modelo de CPW) |
| 3.5.2.2 | CSRR em microstrip | Aba Linha CSRR + gap: eqs. 3.86–3.92, modelo melhorado, geometria inicial, layout |
| 3.5.2.3–3.5.2.4 | Acoplamento entre ressonadores, orientação | Pendente |
| 3.5.2.5 | OSRR e OCSRR em CPW | Extração (G.12–G.19); projeto pendente |
| 3.5.2.6 / Ap. H | Síntese por space mapping (ASM) | Pendente: exige laço com simulador EM |
| 3.5.3 | Linhas híbridas | Pendente |
| Ap. G | Extração de parâmetros | Implementado: 7 modelos, busca entre candidatos, caso balanceado |

## Notas de leitura
- Eq. 3.72: o fator (1 + M²/2LpLs) aparece ao quadrado (confirmado pelas Tabelas 3.1 e 3.2).
- Eq. 3.66: ωCR − ωCL = 2ωR e ωCL·ωCR = ω0², o que permite a síntese direta pela banda.
- CSRR + gap balanceada: o método do Ap. G não determina L (χ(fo) = 0); usa-se ZB em βl = π/2.
