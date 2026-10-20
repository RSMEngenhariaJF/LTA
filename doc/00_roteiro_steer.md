# Microwave and RF Design — Roteiro (Steer)

Fonte: M. Steer, *Microwave and RF Design: A Systems Approach*, SciTech Publishing, 2010.

| Seção | Tema | Situação no LTA (v2.1) |
|---|---|---|
| 5.13 | Linhas diferenciais: impedâncias de modo diferencial e comum | Z0d = 2·Z0o e Z0c = Z0e/2 na análise de linhas acopladas |
| 6.10 | Mudança de impedância de referência | Renormalização de .s2p por porta (aba Caracterização) |
| 6.11 | Calibração TRL | Engen-Hoer com thru de comprimento zero; exemplo sintético (aba Caracterização) |
| 7.3 | Terminações e atenuadores | Pads T e Pi com Z01 ≠ Z02, Kmín, E24/E96, potência (aba Atenuadores) |
| 7.13 | Híbridas e divisores concentrados | Implementação "Concentrada" na aba Híbridas e divisores |
| 8.6 | Redes de três elementos com Q | Redes Pi e T na aba Casamento (Ex. 8.7) |
| 9.10 | Acoplador com capacitores de compensação | Opção na aba Linhas acopladas (comparável à LTA distribuída) |
| 10.14 | Filtro combline | Novo tipo na aba Filtros com linhas, com layout |
| 11.5-11.7 | Estabilidade, ganhos e casamento de amplificadores | Aba Amplificador (pHEMT da Fig. 11-8 embutido) |
| 13.2 | Balun de Marchand miniaturizado | Opção na aba Balun |
| 2-4 | Sistemas, modulação, linhas | Fora do escopo ou já coberto pelo Pozar e pelo Garg-Bahl |
| 11.8-11.10 | Ruído, polarização, amplificadores de potência e distribuídos | Pendente |
| 12 | Osciladores | Pendente |

## Notas
- Tabela 11-7: em 5 GHz o livro traz GTU,max = 16,00 dB; com as S da Fig. 11-8 dá 17,00 dB (os demais pontos conferem).
- Ex. 8.7: X1 = 16,4 Ω no livro; o valor coerente com RV = 4,95 Ω e com Xa = 14,933 Ω é 16,58 Ω.
- Combline: as diferenças de ~0,6% em relação ao livro vêm do protótipo (o LTA usa ε = 0,1 exato, 0,0432 dB).
- Acoplador com capacitores nas pontas: diretividade excelente em f0, mas em banda estreita. A LTA distribuída (8 células) mantém mais de 60 dB numa oitava.
- Marchand miniaturizado: com o acoplamento fixo não há casamento exato (quatro condições para os capacitores); liberando o acoplamento, a solução fica exata com 0,2-1 dB a mais de acoplamento.
