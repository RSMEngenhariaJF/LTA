# Aplicações de linhas metamateriais — Roteiro (capítulo 4)

Fonte: F. Martín, *Artificial Transmission Lines for RF and Microwave Applications*, Wiley, 2015, cap. 4.

> Documento de trabalho. Equações conferidas contra o PDF; diferenças em relação a valores publicados estão anotadas.

| Seção | Tema | Situação no LTA (v1.5) |
|---|---|---|
| 4.2.1 | Componentes de banda larga | Limites de banda (eqs. 4.8–4.15); rat-race com −90° LH (Fig. 4.5) |
| 4.2.2 | Dual-band e multibanda | RH + PLH, CRLH canônica, CSRR + gap (com trechos de linha), E-CRLH quad-band; branch-line, Wilkinson e divisor em Y |
| 4.2.3.1–4.2.3.2 | Filtros rejeita-faixa e supressão de espúrios com SRR/CSRR | Pendente |
| 4.2.3.3 | Filtros e diplexers com células RH/LH alternadas | Pendente |
| 4.2.3.4 | Passa-faixa com células híbridas (CSRR + gap + stub) | Implementado (eqs. 4.38–4.46), periódico |
| 4.2.3.5 | Passa-altas com CRLH balanceada | Coberto pela aba CRLH (banda acima de f_CL) |
| 4.2.3.6 | Banda larga com OSRR/OCSRR (Chebyshev) | Protótipo e inversores; síntese OSRR/OCSRR pendente (falta modelo de CPW) |
| 4.2.3.7 | Passa-baixas elíptico com OCSRR | Pendente |
| 4.2.3 (4.49–4.50) | Dual-band por transformação de frequência | Implementado |
| 4.2.4 | Antenas leaky-wave | Varredura (aba CRLH), nº de células (eq. 4.51), largura de feixe |
| 4.2.5 | Amplificadores distribuídos, filtros ativos | Fora do escopo (circuitos ativos) |
| 4.2.6, 4.3 | Sensores, antenas multibanda, RFID, códigos de barras | Pendente |

## Notas de leitura
- Topologia da E-CRLH (Fig. 3.25c) identificada pelos valores publicados na p. 267 (o texto extraído não a descreve).
- Tabela 4.2: os valores não satisfazem as condições declaradas (fz = 1,80 GHz, Zp(f0) = j33,7 Ω); a solução das equações dá FBW de 9,5% em 1,001 GHz.
- Eq. 4.50: ω0 = √(ω1ω2) e w0 = (ω2 − ω1)/ω0 levam ω1 → −1 e ω2 → +1, que a resposta passa-faixa (simétrica em ±ω) leva ao centro da banda.
