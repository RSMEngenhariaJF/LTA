# Linhas de Transmissão Artificiais Periódicas — Roteiro

Fonte: F. Martín, *Artificial Transmission Lines for RF and Microwave Applications*, Wiley, 2015, cap. 2 (pp. 47–118).

> Documento de trabalho (rascunho). Equações e conclusões devem ser conferidas contra a fonte antes de uso em processo controlado.

## Tópicos do capítulo e potencial de implementação

| Seção | Tema | Equações-chave | O que dá para implementar |
|---|---|---|---|
| 2.2 | Análise de Floquet / harmônicos espaciais | (2.1)–(2.9), (2.11) | Diagrama β×ω com harmônicos espaciais; condição de Bragg f = v_p/(2l); ângulo de radiação leaky-wave cos θ = β_n/k |
| 2.3.1 | Relação de dispersão via ABCD | cosh(γl) = (A+D)/2 (2.21) | Núcleo do pacote: dispersão (α, β) de qualquer célula unitária a partir da matriz ABCD |
| 2.3.2 | Impedância de Bloch | (2.26), (2.27), (2.30)–(2.32) | Z_B± para células simétricas e assimétricas (T vs. L) |
| 2.3.3 | Efeito da assimetria | (2.34)–(2.37), Fig. 2.8 | Reproduzir S11/S21 de 20 células T vs. L (L = 2,5 nH, C = 1 pF) |
| 2.3.5 | Impedância iterativa | (2.40)–(2.42) | Verificação numérica Z_I = Z_B |
| 2.4.1–2.4.2 | Teoria de modos acoplados (CMT) | (2.48)–(2.56) | Solver numérico das equações acopladas a partir de Z₀(z) → S-parâmetros |
| 2.4.3–2.4.4 | Soluções analíticas da CMT | (2.76), (2.84)–(2.89) | f_max, \|S11\|max = tanh(\|k\|L), \|S21\|min = sech(\|k\|L), banda rejeitada |
| 2.4.5 | S11 ≈ transformada de Fourier de K(z) | (2.95)–(2.97) | Estimativa rápida do número/posição de bandas rejeitadas |
| 2.4.6 | Exemplos: K(z) senoidal e quadrado | (2.100)–(2.105) | Reproduzir Figs. 2.11 e 2.14 (RO3010, l = 1 cm, 10 células) |
| 2.5.1 | Aplicações EBG | (2.106)–(2.115) | Refletores, ressonador de Bragg, filtro "wiggly-line", linha de atraso chirpada (CDL) |
| 2.5.2 | Linhas com carga capacitiva (slow-wave) | (2.116)–(2.126) | Dispersão distribuída vs. concentrada; projeto de C, L, C_ls a partir de f_c, Z_B, f_s; reproduzir Fig. 2.29 |
