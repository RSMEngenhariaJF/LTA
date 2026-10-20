# Microstrip Lines and Slotlines — Roteiro (Garg, Bahl e Bianchi)

Fonte: R. Garg, I. Bahl e M. Bianchi, *Microstrip Lines and Slotlines*, 4ª ed., Artech House, 2024.

| Seção | Tema | Situação no LTA (v2.0) |
|---|---|---|
| 1.4.1 | Transição coaxial-microstrip (SMA) | Ls + Cs, predefinidos pela Tabela 2.5, compensação, fixture e de-embedding (aba Conectores) |
| 1.4.2 | Transição guia-microstrip | Pendente |
| 1.5 | Medidas: εr, Z0, εeff, α | Duas linhas, uma linha, dois ressonadores, Q0 → α, VSWR → Z0 (aba Caracterização) |
| 2.4.1, 2.4.5 | Perdas e equações de projeto, espessura do metal | Espessura opcional no substrato (eqs. 2.124-2.127) |
| 2.4.2 | Potência média e de pico | Aba Tolerâncias, potência e orçamento |
| 2.4.3 | Tolerâncias | Sensibilidades, pior caso, ΔW exigido e Monte Carlo |
| 2.4.4 | Anisotropia | Pendente |
| 2.4.6 | Faixa de operação: fT, modo transversal, Qr | Avisos em todas as abas; Qr nos ressonadores |
| 2.4.7 | Modelo RLGC | A partir de S (aba Caracterização) |
| 2.5 | Microstrip suspensa, multicamada, enterrada | Pendente |
| 3.4 | Descontinuidades: ponta, gap, degrau, curva, T, cruzamento, entalhe, via | Ponta, gap, degrau, curva, entalhe e via; T e cruzamento só em gráficos no livro |
| 3.8 | Compensação de degrau, curva e T | Chanfro ótimo nos layouts (eqs. 3.115-3.116) |
| 4.3-4.5 | Linhas acopladas: Kirschning-Jansen, dispersão, espessura, perdas, tolerâncias | Modelo padrão; Wheeler por modo; pior caso de acopladas |
| 4.5.4 | Camada dielétrica sobre as linhas (overlay) | Alternativa à compensação LTA: pendente |
| 4.7 | Multicondutoras | Pendente (Lange usa a aproximação de pares do Pozar) |
| 4.8.2 | Ponta aberta acoplada | Usada no filtro de linhas acopladas e na compensação LTA |
| 4.9 | Balun de Marchand | Aba Balun, com versão LTA |
| 5.3 | Slotline (fórmulas fechadas) | Tipo "Slotline" na aba de linhas convencionais |
| 5.4-5.7 | Descontinuidades, variantes e transições de slotline | Pendente |
| 6 | DGS | Aba DGS: LC, linha periódica, passa-baixas, extração π |
| 7.3-7.4 | CPW/CPS: espessura e perda no condutor | Implementadas |
| 7.7-7.9 | Transições e descontinuidades coplanares | Pendente |
| 8, 9 | Metamateriais CRLH e SIW | Já cobertos pelo livro do Martín |
| 10 | Modelagem numérica | Fora do escopo |

## Notas
- Kirschning-Jansen e Garg-Bahl diferem em até ~6% em Zo para S/h pequeno; o primeiro é o mais preciso (≤ 0,6%).
- A eq. 4.88b desta edição traz 1,5·W/h e tanh(10S/h) no modelo de Garg-Bahl; o LTA usa 2,53·W/h e tanh(8S/h), como na forma do artigo de 1979 reproduzida em outras fontes. A diferença é pequena e o padrão agora é Kirschning-Jansen.
- No filtro LTA, a ponta aberta acoplada (Δle ≠ Δlo) precisa entrar na compensação; sem isso, a espúria em 2f0 volta a −19 dB.
- A transição SMA predefinida reproduz só o VSWR da Tabela 2.5; para projetos reais, extraia Ls e Cs do conector usado.
