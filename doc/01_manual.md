# LTA – Linhas de Transmissão Artificiais Periódicas: Manual

> **Situação:** rascunho técnico, não revisado. Antes de usar os resultados em decisão de projeto ou processo controlado, eles precisam de revisão por pessoa qualificada diferente de quem redigiu e de validação por simulação EM e medição.

Base teórica: F. Martín, *Artificial Transmission Lines for RF and Microwave Applications*, Wiley, 2015, cap. 2. As equações citadas (ex.: eq. 2.21) seguem a numeração do livro.

## 1. Instalação e execução

Requisitos: Python 3.10 ou superior, `numpy` e `matplotlib` (o `tkinter` já vem com o Python para Windows).

```
pip install -r requirements.txt
python main.py            (ou dar duplo clique em executar.bat)
```

Todos os caminhos são relativos à pasta `linhas_periodicas/`. Para levar o projeto a outro computador, basta copiar a pasta inteira.

Testes de validação (comparam com valores publicados no livro):

```
python -m unittest discover testes -v
```

## 2. Navegação e ferramentas

A janela mostra uma ferramenta por vez, escolhida pelo menu:

- **Arquivo:** substrato, pasta de saídas, sair.
- **Ferramentas:** ferramentas agrupadas por capítulo do livro.
- **Exportar:** abre a ferramenta de exportação.
- **Ajuda:** material didático, manual técnico, atalhos e "Sobre".

O substrato ativo fica sempre visível numa barra abaixo do menu. Os atalhos são Ctrl+1…9 para as ferramentas, na ordem do menu, e Ctrl+E para exportar.

**Para acrescentar um capítulo ou ferramenta:** crie a classe (subclasse de `AbaBase`) e inclua-a em `GRUPOS`, em [lta/gui/app.py](../lta/gui/app.py). Um capítulo novo é um grupo novo: menu, atalhos e título são gerados automaticamente.

| Aba | Função | Seções do livro |
|---|---|---|
| **Substrato** | εr, h, t, tan δ, σ e rugosidade usados por todas as abas; calculadora microstrip Z0 ↔ W com εeff(f), Z0(f), δ, αc e αd numa frequência | cap. 1 (eqs. 1.68, 1.79, 1.83, 1.90); Hammerstad-Jensen (1980) |
| **Célula ABCD** | Dispersão (βl, αl), impedância de Bloch e S de N células para topologias T, Π, L. O ramo série é Ls+Cs e o paralelo é Cp‖Lp, o que cobre linhas RH, LH e CRLH | 2.3 (eqs. 2.21, 2.26, 2.40) |
| **Carga capacitiva** | Projeto de linha slow-wave: (fc, ZB, fs) ou (fase em f0, ZB, fs) → L, C, Cls → geometria microstrip com patches; atraso de grupo e atenuação por célula | 2.5.2 (eqs. 2.116–2.126) |
| **EBG / Modos acoplados** | Perfil Z0(z) senoidal ou quadrado → W(z) → S por CMT numérico, cascata ABCD e solução analítica; furos no plano de terra (geometria e Bragg); atraso de grupo | 2.4 e 2.5.1 (eqs. 2.48–2.105) |
| **Floquet / harmônicos espaciais** | Célula EBG (senoidal/quadrada), linha com carga capacitiva ou linha uniforme → βl desdobrado na zona estendida, harmônicos βn, região de onda rápida, ângulo do feixe leaky-wave, vp,n e vg | 2.2 (eqs. 2.6–2.11) |
| **Aplicações EBG (CDL, ressonador)** | Linha de atraso chirpada: (f0, Δf, ψ) → l0, C, L → perfil Z0(z) → S11 e atraso de grupo. Ressonador de Bragg: dois refletores senoidais + cavidade Lm (eq. 2.106) → pico, QL e Q0 | 2.5.1 (eqs. 2.106, 2.109–2.115; 1.74–1.75) |
| **CRLH / D-CRLH** (cap. 3) | Síntese balanceada por banda (f_CL, f_CR, ZB) ou por fases em duas frequências; análise de LR, CR, LL, CL quaisquer. Dispersão com sinal (LH/RH), ZB, S, εeff e μeff, fase e feixe leaky-wave do modo fundamental | 3.4.2–3.4.3 (eqs. 3.51–3.66) |
| **CRLH microstrip (CL-loaded)** (cap. 3) | LR, CR, LL, CL → linha hospedeira, interdigitais em série e stubs com via em derivação; layout com vias | 3.5.1 (Fig. 3.27) |
| **Linha CSRR + gap (microstrip)** (cap. 3) | Projeto balanceado (f0, ZB, fz) do circuito da Fig. 3.38b, modelo melhorado (Fig. 3.40), raio do CSRR e interdigital estimados; layout com o CSRR gravado no terra | 3.5.2.2 (eqs. 3.86–3.92) |
| **Extração de parâmetros** (cap. 3) | Lê .s2p (EM ou medição) e extrai o circuito de CSRR, CSRR + gap, SRR, SRR + fita, OSRR, OCSRR e CRLH CL-loaded; carta de Smith e sobreposição dados × modelo | Apêndice G (eqs. G.1–G.19) |
| **Linhas dual-band e multibanda** (cap. 4) | Fase e impedância em duas (ou quatro) frequências: linha convencional + PLH, CRLH canônica, CSRR + gap com zero de transmissão e trechos de linha, E-CRLH quad-band; verificação nas frequências de projeto | 4.2.2 (eqs. 4.20–4.37) |
| **Híbridas e divisores** (cap. 4) | Branch-line, rat-race, Wilkinson e divisor com inversor, montados como redes de N portas; versão com linhas artificiais contra a convencional; faixas pelos critérios da Tabela 4.1 | 4.2.1–4.2.2 (eqs. 4.2–4.19, Figs. 4.5 e 4.16) |
| **Filtros passa-faixa e dual-band** (cap. 4) | Protótipos Butterworth/Chebyshev; passa-faixa com inversores e ressonadores; transformação dual-band; filtro periódico com células híbridas de CSRR | 4.2.3 (eqs. 4.38–4.50) |
| **Linhas sintonizáveis** (cap. 5) | Varactor, HBV e RF-MEMS (pull-in); linha com SRR e varactor (rejeita-faixa sintonizável); CRLH balanceada sintonizável como defasador (FoM) e antena leaky-wave de varredura em frequência fixa | 5.2–5.3 (eqs. 5.1–5.12) |
| **Linha não linear (NLTL)** (cap. 5) | Sólitons (atraso × amplitude, forma), decomposição de pulsos e multiplicação de frequência por simulação no tempo | 5.4 (eqs. 5.13–5.25) |
| **Linha diferencial (filtro de modo comum)** (cap. 6) | Par diferencial com CSRR/DS-CSRR: zero do modo comum, banda máxima de rejeição (estrutura infinita), filtro de N células para os modos comum e diferencial, banda × ordem | 6.3 (Fig. 6.12, eq. 6.11, Tabela 6.1; eqs. 3.93–3.97) |
| **Linhas de atraso MIW/EIW** (cap. 6) | Cadeias de SRR/CSRR acoplados: dispersão, banda, fator de onda lenta; etiqueta RFID por reflectometria (bits, ressonadores) e ecos de um pulso | 6.2 (eqs. 6.1–6.10) |
| **Guia integrado ao substrato (SIW)** (cap. 6) | Largura, vias e banda monomodo pelo corte do TE10; dispersão, λg, Z_TE, perda dielétrica; layout com vias e transições em cunha | 6.5 (eqs. 6.31–6.34) |
| **Linha com capacitores negativos** (cap. 6) | Linha ENZ não Foster: C' efetiva, estabilidade, Z, ε equivalente, vp e vg superluminais | 6.4.2 (Fig. 6.32) |
| **Exportar** | Layout para KiCad, DXF, Gerber, GDSII, script HFSS e macro CST; parâmetros S em Touchstone e CSV | — |

Cada aba de projeto tem botões **Exemplo** que reproduzem figuras do livro, **Enviar p/ exportação** (manda o layout para a aba Exportar) e **Exportar S2P**.

### Banda larga: redes passa-tudo e banda útil

- **Topologias Treliça (X) e T em ponte** (aba Célula ABCD), baseadas em Martín (2015), seção 6.4, pp. 429–439:
  - Na treliça, Z_B = √(Zs·Zp) é constante quando Zs·Zp é real (ex.: Ls nos braços e Cp nas diagonais). O resultado é uma rede passa-tudo, com |S21| = 0 dB e S11 ≈ 0 em qualquer frequência. É uma rede balanceada (diferencial).
  - No T em ponte (Z3 = 4·Z2, montada automaticamente), Z_B = √(2·Z1·Z2). Tem terminação simples (single-ended).
  - Os botões *Ex. passa-tudo (treliça)* e *Ex. T em ponte* carregam os casos das Figs. 6.25a e 6.31 do livro.
- **Banda útil:** as abas Célula, Carga capacitiva e EBG informam as faixas em que |S11| ≤ −15 dB e |S21| ≥ −1 dB. Os critérios são editáveis na aba Célula.
- **Regra prática para escadas LC:** S11 ≤ −15 dB, para qualquer número de células, é garantido até f ≈ 0,55·fc (fc ≥ 1,8·fmax).
- **EBG:** multiplicar K(z) por `janela(z)` ou `hann(z)` no perfil personalizado reduz a ondulação e amplia a faixa útil abaixo da banda rejeitada. No exemplo do livro, ela passa de 2,36 para 3,68 GHz.

### Floquet e radiação leaky-wave

Aba **Floquet / harmônicos espaciais** ([floquet.py](../lta/core/floquet.py)).

- O βl da eq. 2.21 está na zona reduzida [0, π]. A ferramenta o **desdobra** na zona estendida: escolhe, entre 2πk ± βl, o valor contínuo e não decrescente. Isso funciona com perdas e com bandas rejeitadas mais estreitas que o passo de frequência. Vale para estruturas de mão direita; linhas LH/CRLH (cap. 3) precisarão de outro critério.
- **Gráficos:**
  - Diagrama de Floquet em zona estendida (Fig. 2.2), com βn·l = βl + 2πn (eq. 2.6), os ramos simétricos (pontilhados) e a região de onda rápida |βn| < k0 sombreada. Os trechos que radiam ficam destacados.
  - Atenuação por célula e βl estendido.
  - Ângulo do feixe pela eq. 2.11 (cos θ = βn/k0, θ a partir do eixo da linha; 90° = broadside).
  - vp,n/c (eq. 2.8) e vg/c (eq. 2.7).
- **Texto:** bandas rejeitadas, Bragg (eq. 2.9), estimativas de início de radiação e broadside do harmônico n = −1, e faixas de radiação de cada harmônico (para trás, βn < 0, espectral; para a frente, βn > 0, não espectral; eq. 2.13).
- A radiação **não** entra em α. A ferramenta indica onde ela ocorre, mas não quanto se perde.
- **Exemplos:** EBG da Fig. 2.11; linha uniforme (Fig. 2.1: retas deslocadas de 2π/l); linha com carga capacitiva.

### Aplicações EBG: CDL e ressonador de Bragg

Aba **Aplicações EBG** ([aplicacoes_ebg.py](../lta/core/aplicacoes_ebg.py)).

- **Linha de atraso chirpada (CDL), seção 2.5.1.5:**
  - **Entradas:** f0, Δf e ψ (inclinação do atraso de grupo, em ns/GHz). O sinal de ψ define o do chirp: ψ < 0 reflete primeiro as frequências altas, que são as mais atenuadas.
  - **Projeto:** l0, C e L pelas eqs. 2.115. L é arredondado para N·l0, para que as portas fiquem em Zref, e C é recalculado pela eq. 2.112 para manter Δf. Também se pode fixar N, como no livro (N = 50).
  - **Perfil:** Z0(z) = Zref·exp(A·W(z)·sin(2πz/l0 + Cz² − CL²/4)) (eq. 2.109). Janela gaussiana assimétrica do livro, W = exp(−4((z − L/4)/L)²), Hamming ou nenhuma.
  - **Resultados:** perfil, frequência refletida em cada posição (eq. 2.111), |S11| e |S21|, atraso de grupo de S11 com a reta de projeto e ajuste linear do ψ simulado.
  - **Exemplo do livro (Fig. 2.26):** l0 = 6,39 mm (livro 6,4), C = −2052 m⁻² (livro −2080), ψ simulado sem perdas −0,48 ns/GHz (projeto −0,46). Com perdas, a perda de reflexão média é ≈ 3,4 dB (livro: ≈ 3 dB).
  - Com dispersão ligada, a banda real se comprime e ψ muda (≈ +19% no exemplo), porque as eqs. 2.111–2.115 supõem εeff fixo.
- **Ressonador de Bragg, seção 2.5.1.2:**
  - **Estrutura:** dois refletores senoidais de largura modulada, sintonizados em f0, separados por uma cavidade uniforme de comprimento Lm.
  - **Projeto:** as fases φref1 (S22 do refletor 1) e φref2 (S11 do refletor 2) vêm do CMT em f0. Lm vem da eq. 2.106, φ1 + φ2 − 2βLm = 2πn, com o n que dá o menor Lm positivo.
  - **Resultados:** pico de |S21|, perda de inserção, QL = f0/BW−3dB (eq. 1.75) e Q0 = QL/(1 − 10^(−IL/20)) (eq. 1.74).
  - **Exemplo (Fig. 2.21, 4,2 GHz em RO3010):** pico exatamente em f0 e QL = 132 com perdas. O livro mede Q = 129,2, mas com refletores de furos no terra, que não têm modelo de circuito aqui.
- Os dois casos exportam o layout (largura modulada) e o Touchstone.
- **Filtro wiggly-line (seção 2.5.1.3):** não implementado. Ele exige um modelo de linhas acopladas, que o programa ainda não tem.

### Escala linear e Fig. 2.8 (aba Célula ABCD)

- O campo **Escala de |S|** alterna entre dB e linear.
- O botão **Exemplo Fig. 2.8** reproduz a figura do livro: 20 células T contra L, L = 2,5 nH, C = 1 pF, em escala linear.
- Quando se comparam células T e L de linha convencional, o texto informa ω0 = 2/√(LC) e o limite de equivalência ω ≪ ω0/2 (eqs. 2.35–2.37).

### Perdas e dispersão (modelo físico da linha)

As abas **EBG / Modos acoplados** e **Carga capacitiva** têm a caixa **Modelo físico da linha**, com duas opções independentes. Desligadas, o cálculo é o quase-estático sem perdas das versões anteriores, com resultado idêntico.

- **Perdas:** α = αc + αd, e a constante de propagação passa a ser γ = α + jβ ([microstrip.py](../lta/core/microstrip.py)).
  - **Dielétrico:** eq. 1.68.
  - **Condutor:** eq. 1.90 para W ≥ 2h e Pucel/Gupta (1979) para tiras mais estreitas. Rs = 1/(σδ) (eqs. 1.79 e 1.83), com correção de espessura finita (t comparável a δ) e de rugosidade (Hammerstad: K = 1 + (2/π)·atan(1,4·(Δ/δ)²)).
  - σ, t, rugosidade e tan δ vêm da aba Substrato.
- **Dispersão:** εeff(f) por Kirschning-Jansen (1982), válido para εr ≤ 20 e f·h ≤ 25 GHz·mm; Z0(f) por Bianco et al. (1978).
  - Na EBG do exemplo Fig. 2.11, o centro da banda passa de 5,82 para 5,64 GHz. O livro (solver EM) dá 5,66 GHz.
  - Com dispersão, o Z0 das pontas da linha sobe com a frequência (≈ 54 Ω a 10 GHz em RO3010). A cascata ABCD, referida a Zref fixo, mostra esse descasamento; o CMT usa como referência o Z0 local das pontas. Por isso, as curvas de |S11| se afastam em alta frequência.
- **EBG:** o CMT numérico (γ complexo nas equações acopladas) e a cascata ABCD aceitam Z0, εeff e α por frequência. O texto informa a **perda dissipativa** −10·log(|S11|²+|S21|²) na banda passante e os centros de banda do cálculo numérico. O último gráfico mostra α(f) e εeff(f) da linha de referência.
- **Carga capacitiva:**
  - A linha de alta impedância recebe α(f) e εeff(f) da sua largura (`celula_microstrip`).
  - A perda do patch é desprezada, e Cls continua concentrada.
  - Gráficos novos: atraso de grupo e atenuação por célula (αl de Bloch).
- **Atraso de grupo:** τg = −dφ(S21)/dω ([abcd.py](../lta/core/abcd.py)). Use pontos suficientes para que a fase não salte mais de π entre amostras.
- As soluções analíticas da seção 2.4 (eqs. 2.76–2.89) continuam sem perdas e sem dispersão.

### Perfil EBG personalizado (função de z)

Na aba **EBG / Modos acoplados**, a perturbação **Personalizado (função)** aceita uma expressão em z (mm, de 0 a L) que define **K(z) [1/cm]**, **Z0(z) [Ω]** ou **W(z) [mm]** ([perfil_personalizado.py](../lta/core/perfil_personalizado.py)).

- **Várias frequências:** cada termo `A*cos(2*pi*z/bragg(f))` em K(z) cria uma banda em f GHz, com |S21|min ≈ sech(A·L/2) (L em cm). Exemplo: `0.4*cos(2*pi*z/bragg(3)) + 0.4*cos(2*pi*z/bragg(5))`.
- **Funções auxiliares:** `bragg(f)` (período em mm para a banda em f GHz), `janela(z)` (Hamming), `hann(z)`, `ret(z,a,b)`, `degrau(z,a)`, além de `sin`, `cos`, `exp`, `log`, `sqrt`, `tanh`, `sign`, `where`, `clip`, `pi`. Variáveis: `z`, `L`, `Zref`.
- No modo K(z), o valor médio de K é removido para a linha terminar em Zref.
- **Modelos prontos** (botão *Carregar*): duas bandas, três bandas com janela, supressão de 2f0 e 3f0, chirp de banda larga, perfil do livro em Z0(z) e largura quase quadrada.
- **Resultados:** CMT, cascata ABCD, estimativa de Fourier de K(z) sobre o |S11| (Eq. 2.95) e lista automática das bandas abaixo de −10 dB.
- **Segurança:** a expressão é verificada pela árvore sintática antes de ser avaliada. Só números, operadores, comparações e as funções listadas são aceitos.

### Pré-visualização da estrutura

O lado direito das abas de projeto tem duas sub-abas: **Estrutura (pré-visualização)** e **Resultados**.

- A prévia **atualiza enquanto você edita** os campos, cerca de 0,4 s depois de parar de digitar, sem precisar clicar em Calcular. Ao clicar em Calcular, a tela passa para Resultados.
- **Vistas:**
  - placa inteira, com comprimento e largura totais e a célula do zoom destacada;
  - zoom da célula unitária, com cotas marcadas por letras e os valores no quadro do canto;
  - corte transversal A, com substrato, terra (e abertura do furo, se houver) e trilha. A espessura do cobre está exagerada para ficar visível.
- **Comparação:** a versão anterior aparece tracejada em cinza, e a caixa de texto lista o que mudou (ex.: `l_patch_mm: 2.77 → 4.35`). O botão **Fixar referência** congela uma versão para comparar várias alterações contra ela; **Soltar referência** volta a comparar com o estado anterior. Estruturas de tipos diferentes não são sobrepostas.
- **Avisos:** especificação inviável (ex.: patch maior que o período) aparece em vermelho sobre a última estrutura válida. Trilhas com menos de ~0,1 mm geram alerta de fabricabilidade.
- Na aba **Célula ABCD**, a prévia é o **esquemático** da célula (T, Π ou L) com os valores de cada elemento.
- A aba **Exportar** usa a mesma prévia para a estrutura selecionada.

### Capítulo 3: linhas metamateriais

As quatro ferramentas do grupo **Cap. 3 – Metamateriais** seguem o fluxo do livro: **circuito → geometria inicial → simulação EM ou medição → extração → correção**. Isso é necessário porque, no capítulo 3, a ligação entre circuito e geometria não tem fórmula fechada confiável; o livro a fecha com EM e extração de parâmetros (Apêndice G) ou com *space mapping* (Apêndice H).

- **CRLH / D-CRLH** ([crlh.py](../lta/core/crlh.py)):
  - **Síntese balanceada por banda:** ωCR − ωCL = 2ωR e ωCL·ωCR = ω0² (eq. 3.66), com ZB = √(LR/CR) = √(LL/CL) (eq. 3.63).
  - **Síntese por fases:** βl = ω/ωR − ωL/ω (eq. 3.59 balanceada), para linhas dual-band. Vale com |βl| pequeno por célula.
  - **Sinal de β:** negativo onde a reatância série é negativa. É a banda LH (eq. 3.60 generalizada). O feixe leaky-wave do modo n = 0 varre de backfire (LH) a endfire (RH), passando por broadside na transição.
  - **Exemplo não balanceado:** reproduz a Fig. 3.21, com gap de ωs a ωp e ZB nula em ωs e com polo em ωp.
- **CRLH microstrip CL-loaded** ([cl_loaded.py](../lta/core/cl_loaded.py)):
  - A linha hospedeira de Z = √(LR/CR) e período p dá LR.
  - O capacitor interdigital (fórmula de Bahl) dá CL. Os capacitores das pontas valem 2CL.
  - O stub em curto (Z0 escolhida) mais a via (Goldfarb e Pucel) dão LL.
  - Os avisos indicam interdigital grande demais para o período ou stub longo. Isso é a limitação real de banda larga em uma camada, que o livro discute na seção 3.5.1.
  - O exemplo LR = 2,45 nH, CR = 0,5 pF, LL = 3,38 nH, CL = 0,68 pF dá f0 = 3,9 GHz, como a linha da Fig. 3.28.
- **Linha CSRR + gap** ([csrr.py](../lta/core/csrr.py)):
  - **Projeto do circuito:** Cc = L/ZB², Lc = 1/(ω0²Cc), Cg = 1/(ω0²L) e C = 1/(ωz²Lc) − Cc (eqs. 3.86, 3.91, 3.92). L vem da linha no período escolhido.
  - **Modelo melhorado (eqs. 3.87–3.88 invertidas):** separa a capacitância série física do gap, Cs = Cg²/(C + 2Cg), de Cpar = Cg·C/(C + 2Cg). Na Tabela 3.3, Cs = 0,014 pF e Cpar = 0,28 pF: o Cg do modelo é quase todo capacitância paralela, como o livro argumenta.
  - **Raio do CSRR:** estimativa quase-estática da seção 3.3.2.1 (Lc = Lo/4 por CPW; Cc por dualidade com o SRR). Contra o livro, fo sai 16% baixa na Fig. 3.39.
  - **Fator de calibração:** fo real / fo estimada, obtido de uma simulação ou medição pela aba Extração. Com a calibração, o raio da Fig. 3.39 (3,661 mm) é reproduzido.
- **Extração de parâmetros** ([extracao.py](../lta/core/extracao.py)):
  - **Condições do Apêndice G:** S11 sobre os círculos r = 1 (ramo em derivação abre) e g = 1 (ramo série se anula), zero de transmissão, zero de reflexão e βl = ±π/2.
  - **Candidatos:** os cruzamentos com os círculos e os zeros do equivalente T ou π da célula medida (que pegam pontos de tangência). As combinações são testadas e fica a de menor erro rms.
  - **Linha balanceada (fs ≈ fo):** a reatância em fo se anula e não determina L. Usa-se ZB(βl = π/2) = |X do braço série|, que é exato para célula simétrica sem perdas.
  - **Validação:** o modelo SRR + fita reproduz as transformações das eqs. 3.71–3.74 entre as Tabelas 3.2 e 3.1. Na eq. 3.72, o fator (1 + M²/2LpLs) aparece ao quadrado; o texto extraído do PDF não mostra o expoente, mas só assim a Tabela 3.1 é reproduzida.
  - **Fontes de dados:** arquivo .s2p (Touchstone v1: RI, MA ou DB, em qualquer unidade), resultados das outras abas (ex.: `celula_csrr`, `celula_crlh_cl`) ou exemplos com os valores das tabelas do livro.
  - **Erro rms:** < 0,02 indica bom ajuste; > 0,1 indica modelo inadequado ou faixa grande demais.

### Capítulo 4: aplicações

As ferramentas do grupo **Cap. 4 – Aplicações** trabalham no nível de circuito, com o mesmo fluxo do capítulo 3 para chegar ao layout (abas CRLH CL-loaded e CSRR, depois EM e extração).

- **Linhas dual-band e multibanda** ([multibanda.py](../lta/core/multibanda.py)):
  - **Condição geral:** a célula T com fase φ e impedância Z tem braços Xa = Z·tan(φ/2) e derivação Xp = −Z/sen φ (inversão das eqs. 2.30 e 2.33). Com φ = ∓90°, são as eqs. 4.28.
  - **Linha convencional + PLH:** eqs. 4.21–4.23, com a condição 4.24. As fórmulas supõem a PLH homogênea: com uma única célula a fase erra 23° em f1; com 4 células, 2°. A tabela de verificação mostra isso.
  - **CRLH canônica:** solução linear (eqs. 4.29); com Z1 = Z2, a célula é balanceada em √(ω1ω2) (eqs. 4.31–4.32).
  - **CSRR + gap:** o ramo série sai direto (eqs. 4.33a–b). O ramo em derivação (C, Lc, Cc) usa as fases e o zero de transmissão (eqs. 4.34 generalizadas). Com trechos de linha nas pontas, φcélula = φ − 2φL (p. 259). Os dois exemplos do livro são reproduzidos a 3 algarismos.
  - **Quad-band E-CRLH:** oito condições (eqs. 4.37), resolvidas por ramo, com 3 incógnitas lineares e 1 por busca de raiz. A topologia da Fig. 3.25c foi identificada pelos valores publicados: braço Lhs/2 + 2Chs + (2Lhp ‖ Chp/2); derivação (Lvp ‖ Cvp) ‖ (Lvs + Cvs). A solução satisfaz exatamente ±90° e 35,35 Ω nas quatro frequências; os valores do livro ficam a até 3,4% dela e erram até 4 Ω.
- **Híbridas e divisores** ([dispositivos.py](../lta/core/dispositivos.py)):
  - **Rede de N portas:** cada linha entra pela sua matriz ABCD na matriz Y nodal; os nós internos são eliminados (redução de Kron).
  - **Dispositivos:**
    - branch-line: entrada 1, direta 2, acoplada 3, isolada 4;
    - rat-race: porta 4 = Σ, porta 1 = Δ, como na Fig. 4.7;
    - Wilkinson;
    - divisor em Y com inversor (Figs. 4.11 e 4.15).
  - **Linhas:** convencional, CRLH dual-band (−90° em f1, +90° em f2) ou PLH de −90° no lugar do ramo de 270° do rat-race (banda larga, Figs. 4.1 e 4.5).
  - **Métricas:** faixa contínua em torno de cada frequência de projeto em que valem |S11| < −10 dB, saídas −3 ± 1 dB, balanço de fase ± 5° e isolamento < −15 dB.
  - **Limites de banda da seção 4.2.1:** dφ/dω da linha dual (eqs. 4.8–4.11) e a fase-limite de uma célula, 0,708π (livro: 0,7π).
- **Antena leaky-wave** (aba CRLH): com α/k0 obtido por EM, a aba calcula o nº de células para irradiar a fração escolhida da potência (eq. 4.51) e a largura de feixe pela aproximação de abertura (Jackson e Oliner).
- **Filtros** ([filtros.py](../lta/core/filtros.py)):
  - **Protótipos:** Butterworth e Chebyshev (inclusive ordem par, com o último inversor K = Z0·√g(n+1)).
  - **Passa-faixa:** inversores ideais e tanques Ceq = g/(FBW·ω0·Z0) (eqs. 4.40–4.42).
  - **Dual-band:** transformação das eqs. 4.49–4.50, que leva cada tanque a (L' ‖ C') ‖ (L'' + C''), a derivação da E-CRLH.
  - **Filtro periódico com células híbridas de CSRR:**
    - Cg pela eq. 4.46; Ld, C, Lc e Cc pelas eqs. 3.86 e 4.43–4.45;
    - g de elementos iguais com 3 dB em Ω = 1 dá 1,521 para n = 3, como no livro;
    - a solução atende às condições declaradas e dá FBW de 9,5% centrada em 1,001 GHz;
    - os valores da Tabela 4.2 dão fz = 1,80 GHz e Zp(f0) = j33,7 Ω. Provavelmente incluem o ajuste de layout descrito no livro, então a diferença em relação a eles é esperada.

### Capítulo 5: linhas sintonizáveis e não lineares

- **Dispositivos de sintonia** ([sintonia.py](../lta/core/sintonia.py)):
  - **Varactor:** Cj0/(1 + V/φ)^γ, com γ entre 1/3 e 1/2, conforme a nota 1 do livro.
  - **HBV:** CM·sech(kV), capacitância par.
  - **RF-MEMS de placas paralelas:** C = ε0A/g; a tensão V(g) pela eq. 5.4; o pull-in Vp = √(8k·g0³/(27ε0A)) em g = 2g0/3 (eq. 5.5).
- **Linha com SRR e varactor (VLSRR):**
  - Ceq = CL em série com Cvar (CR é desprezada, como no texto).
  - O tanque refletido na linha (eqs. 5.6) tem Ls' = 4M²Ceqω0² = 4M²/Lr, Cs' = Lr/(4M²ω0²) e R = 4M²ω0²/Rs, com Q = ω0Lr/Rs.
  - Com os valores da Fig. 5.7, sintoniza de 2,64 a 3,80 GHz, com rejeição mais funda no topo da faixa, como o livro descreve.
- **CRLH sintonizável:**
  - CL e LL são multiplicados pelo mesmo fator x (eqs. 5.9–5.11). Na linha balanceada, ZB fica constante e a transição vai para f0/√x.
  - **Defasador:** mostra a defasagem em f_op e a FoM da eq. 5.12, com Q finito nos elementos sintonizáveis.
  - **Antena leaky-wave:** ângulo do feixe em frequência fixa (seção 5.3.2.2).
- **Linha não linear** ([nltl.py](../lta/core/nltl.py)):
  - **Sólitons:**
    - F(V) e G(V) pelas eqs. 5.21 e 5.25;
    - faixa de atrasos possíveis: √(L·C(∞)) a √(L·C(0));
    - atraso × amplitude pela eq. 5.24;
    - forma do pulso por quadratura da eq. 5.20: (T⁴/24L)(dV/ds)² = G(V) − T²V²/(2L).
  - **Simulação no tempo:** integra as eqs. 5.13–5.15 (RK4) com fonte e carga resistivas, vetorizada para simular várias linhas de comprimentos diferentes de uma vez.
  - **Validação interna:** o atraso dos sólitons simulados confere com a eq. 5.24 com erro ≤ 1,5%. O erro cresce para sólitons curtos, como esperado da aproximação contínua (eq. 5.17).
  - **Multiplicador:** com os parâmetros da Fig. 5.44, o máximo do 3º harmônico cai com CM, como no livro: 33 seções para CM = 0,12 pF, contra 27 no livro.
    - A posição exata depende das terminações, que o livro não informa; o padrão aqui é √(L/Co) = 50 Ω.
    - As eficiências não incluem perdas.

### Capítulo 6: outras linhas avançadas

- **Linha diferencial com supressão do modo comum** ([diferencial.py](../lta/core/diferencial.py)):
  - **Modo comum (parede magnética):** cada linha vê Le e o acoplamento Ce ao meio CSRR (2Lc, Cc/2); os CSRRs vizinhos se acoplam por CR (Fig. 3.44b).
  - **Zero de transmissão (eq. 6.11):** fz = 1/(2π√(2Lc(Ce + Cc/2))), que confirma a leitura em meio circuito dos valores da Tabela 6.1 (1,495 GHz para a linha de 1,5 GHz).
  - **Banda máxima de rejeição:** faixa em que nenhum dos dois modos de Bloch da célula de 4 portas propaga (eqs. 3.96–3.97). Reproduz a Tabela 6.1 com no máximo 2,0 pontos de diferença.
  - **Filtro de N células:** simulado como rede; a banda a −20 dB converge para a máxima, e seis células já ficam perto dela, como diz o livro.
  - **Modo diferencial (parede elétrica):** o CSRR não é excitado e a linha (Lo, Co pelo período e pela εeff do modo ímpar) passa com perda desprezível.
- **Linhas de atraso MIW/EIW** ([miw.py](../lta/core/miw.py)):
  - **Dispersão:** ω = ω0/√(1 + 2κ·cos βl), com κ = M/L (MIW; M < 0 no arranjo coplanar, onda regressiva) ou κ = −CM/Cc (EIW, eq. 6.4). Com |M|/L = CM/Cc as duas bandas coincidem, como prevê a dualidade.
  - **Velocidade de grupo e atraso:** vg = l·κ·ω0 no centro (eq. 6.8) e atraso de ida e volta 2/(ω0|κ|) por ressonador (eq. 6.9), conferido contra o atraso de grupo da cadeia simulada.
  - **Etiqueta RFID por reflectometria:** bits = Δτ/τ (eqs. 6.7 e 6.10). Os refletores ficam ao fim de cada segmento de bit, e a aba mostra o envelope do eco de um pulso gaussiano.
  - **Limitação:** pulsos de banda larga em relação à banda MIW chegam atrasados e geram ecos espúrios. É o compromisso entre bits, banda, comprimento e perdas descrito no livro (p. 429).
- **SIW** ([siw.py](../lta/core/siw.py)):
  - **Projeto:** a = c/(2fc√εr) (eq. 6.34), W pela eq. 6.31 (concorda com a fórmula de Deslandes e Wu a 0,6%), d máx. pela eq. 6.32 e s ≤ 2d pela eq. 6.33.
  - **Modo TE10:** β, λg, Z_TE e αd; trecho de guia como linha equivalente.
  - **Layout:** fileiras de vias e transições em cunha (dimensões iniciais, a otimizar em EM).
  - **Fora do escopo:** os filtros e linhas CRLH em SIW (seções 6.5.1–6.5.2) dependem de simulação EM.
- **Linha com capacitores negativos** ([naofoster.py](../lta/core/naofoster.py)):
  - **Modelo:** C'_ef = C' − Cn/d (estável se > 0), Z = √(L'/C'_ef) e ε = c²L'C'_ef < 1, com vp = vg > c (conferido pelo atraso de grupo da rede).
  - **Divergência com o livro:** com os dados da Fig. 6.32, o modelo dá 48,3 Ω e ε = 0,43; o livro cita 50 Ω e 0,3.

## 3. Exportação

### 3.1 Layout

| Arquivo | Ferramenta | Como abrir |
|---|---|---|
| `.kicad_pcb` | KiCad 6/7/8 | Abrir no editor de PCB e pressionar **B** para preencher a zona de terra |
| `.kicad_pcb` | Altium | *File › Import Wizard › KiCad* (versões recentes; confirmar na versão instalada) |
| `.dxf` (R12) | Altium, ADS, HFSS, CST | Altium: *File › Import › DXF/DWG* no PCB, mapeando as camadas |
| `_topo.gtl`, `_terra.gbl`, `_contorno.gko` | Altium (CAMtastic), fabricação | Importar como Gerber RS-274X, em mm |
| `_furos.drl` | Fabricação, Altium | Furação Excellon (mm) das vias metalizadas; só existe se houver vias |
| `.gds` | ADS, HFSS, CST | ADS: *Layout › File › Import › GDSII*. Unidade de banco de dados 1 nm |
| `_hfss.py` | Ansys HFSS / AEDT | *Tools › Run Script*. Cria o projeto com variáveis `h_sub` e `h_ar`, portas de onda, radiação, setup e varredura |
| `_cst.bas` | CST Studio Suite | *Home › Macros* (editar/executar macro VBA). Parâmetros `h`, `t`, `er` e `tand` ficam editáveis |
| `_LEIA-ME.txt` | — | Substrato, camadas e parâmetros de projeto da estrutura |

Camadas (DXF/GDS): `TOPO_COBRE` (GDS 1), `TERRA` (GDS 2), `TERRA_FUROS` (GDS 3, furos e fendas a subtrair do terra, como os anéis do CSRR), `VIAS` (GDS 4) e `PLACA` (GDS 10). No KiCad, as fendas viram áreas de exclusão do preenchimento do terra e as vias viram objetos *via*. No HFSS e no CST, as fendas são subtraídas do terra e as vias são cilindros PEC do terra ao topo.

**Por que não "PCB ASCII" do Altium:** o formato ASCII do `.PcbDoc` não tem especificação pública, e um arquivo gerado por engenharia reversa pode falhar silenciosamente entre versões. As rotas DXF, Gerber e KiCad são documentadas e o Altium as importa.

### 3.2 Parâmetros S

Na aba Exportar, a seção "Parâmetros S calculados" lista o último resultado de cada aba:

- **`.s2p`** (Touchstone v1, formato RI, frequência em GHz, Z0 de referência no cabeçalho). No ADS, use o componente *Data Items › S2P* apontando para o arquivo, que pode ser cascateado com outros blocos no esquemático.
- **`_S.csv`** com |S| em dB e fase em graus de S11, S21, S12 e S22, separados por `;`.

Os parâmetros S vêm de **modelos de circuito ideais**: sem perdas, sem descontinuidades e sem radiação. Servem para estudo de sistema e comparação, mas não substituem o S simulado em EM ou medido.

## 4. Validação (testes automatizados)

| Item | Livro | Ferramenta |
|---|---|---|
| W de 50 Ω em RO3010 (h = 1,27 mm) | 1,2 mm | 1,186 mm |
| Z0 máx/mín, perfil senoidal \|K1\| = 0,64 cm⁻¹ | 75 / 33,3 Ω | 75,1 / 33,3 Ω |
| \|S21\|mín senoidal (eq. 2.86) | −49,6 dB | −49,6 dB (analítico e numérico) |
| Banda rejeitada senoidal (eq. 2.89) | 2,55 GHz | 2,61 GHz |
| \|S21\|mín quadrada, n = 1, 3, 5 | −49,3 / −12,5 / −5,7 dB | −49,3 / −12,5 / −5,7 dB |
| f_max (eq. 2.84) | 5,66 GHz | 5,74 GHz (\*) |
| Centro da banda, Fig. 2.11, CMT com dispersão | 5,66 GHz | 5,64 GHz |
| Profundidade pelicular do cobre, 1 e 10 GHz (nota 24, cap. 1) | 2,06 / 0,66 µm | 2,09 / 0,66 µm |
| CDL (p. 101): l0, C com L = 50·l0 | 6,4 mm, −2080 m⁻² | 6,39 mm, −2052 m⁻² |
| CDL: inclinação do atraso de grupo (CMT sem perdas) | −0,5 ns/GHz (−0,46 com L = 50·l0) | −0,48 ns/GHz |
| Ressonador (p. 93): βLm para φref = −95,84° | 84,15°, n = −1 | 84,16°, n = −1 |
| Ressonador projetado: pico de \|S21\| | f0 | f0 (erro < 0,1%) |
| Floquet, linha uniforme: início da radiação (n = −1) e broadside | c/(l(√εeff+1)), c/(l√εeff) | idem (θ = 90° ± 1°) |
| CRLH balanceada: bordas (eq. 3.66) e ZB(f0) (eq. 3.63) | f_CL, f_CR especificadas; ZB | exatas; 50,00 Ω |
| CRLH/D-CRLH: cos βl analítico (eqs. 3.51, 3.64) contra ABCD | — | erro < 10⁻⁸ |
| Transformação SRR + fita (eqs. 3.71–3.74), Tabela 3.2 → 3.1 | L' 4,91 / Lp' 1,38 / Ls' 1,68 / Cs' 3,90 | 4,95 / 1,37 / 1,70 / 3,83 |
| Idem, sem fita | L' 1,91 / Ls' 0,24 / Cs' 23,01 | 1,91 / 0,24 / 22,76 |
| Extração (7 modelos), ida e volta com os valores das Tabelas 3.1 e 3.3 | valores das tabelas | desvio < 0,5% |
| Extração com perdas no ressonador (R = 5 kΩ) | — | desvio < 2% |
| Linha CSRR, ZB pela eq. 3.90 contra ABCD; eqs. 3.87–3.88 | — | exatas |
| CSRR da Fig. 3.39: fo estimada pela geometria | 2,34 GHz (extraída) | 1,97 GHz (−16%) |
| CRLH CL-loaded: f0 com os elementos do exemplo | 3,9 GHz (Fig. 3.28) | 3,885 GHz |
| Dual-band CSRR + gap, p. 258 (L, C, Cg, Lc, Cc) | 12,5 / 24,9 / 1,25 / 3,38 / 5,10 | 12,50 / 24,91 / 1,251 / 3,376 / 5,096 |
| Idem com trechos de linha, φL(f1) = 15° (p. 259) | 9,45 / 17,9 / 1,01 / 4,85 / 2,95 | 9,45 / 17,92 / 1,005 / 4,852 / 2,956 |
| Inversor CRLH 2,4/3,75 GHz (p. 255): CL, LL | 0,33 pF, 0,83 nH (após parasitas) | 0,338 pF, 0,844 nH |
| Quad-band E-CRLH (p. 267), 8 elementos | valores publicados | a 3,4% (e ±90°, 35,35 Ω exatos) |
| Fase-limite de uma célula dual (p. 239) | 0,7π | 0,708π |
| Rat-race com −90° LH: banda de balanço de fase ±5° | maior que a convencional | 15,6% → 102% (circuito ideal) |
| g periódico de 3 elementos (p. 283) | 1,521 | 1,5214 |
| Chebyshev 0,5 dB n = 3 / 0,1 dB n = 4 (tabelas) | 1,5963 1,0967 / 1,1088 1,3061 1,7703 0,8180 1,3554 | idem (±2·10⁻⁴) |
| Filtro híbrido periódico, f0 = 1 GHz, FBW 10% | 8% medido | 9,5% (circuito) |
| RF-MEMS: pull-in em g = 2g0/3 (eq. 5.5) | 2g0/3 | 0,6667·g0 |
| VLSRR (Fig. 5.7): faixa de sintonia, Q = ω0Lr/Rs | rejeição mais funda no topo | 2,64–3,80 GHz; −32 → −48 dB |
| CRLH balanceada sintonizada: ZB (eq. 5.11) | constante | 50,00 Ω em x = 0,6 a 1,7 |
| NLTL (Fig. 5.40): faixa de atraso dos sólitons | — | 2,174 a 4,216 ps |
| Sóliton simulado no tempo × eq. 5.24 | — | erro ≤ 1,5% |
| Multiplicador NLTL (Fig. 5.44): seções do máximo de 3º harmônico, CM = 0,12 pF | 27 | 33 (cai com CM, como no livro) |
| Linha diferencial: fz (eq. 6.11) com os elementos da Tabela 6.1 | 1,5 / 3 / 6 GHz | 1,495 / 2,906 / 5,931 (CSRR) |
| Linha diferencial: banda máxima de rejeição (Tabela 6.1) | 30,7 / 26,4 / 20,7 / 32,6 / 29,4 / 23,6% | 32,7 / 26,9 / 21,6 / 32,4 / 28,9 / 21,9% |
| MIW: banda com \|M\|/L = 0,0256 (p. 429) | 5,13% | 5,13% (circuito = eq. 6.1 a 10⁻¹⁶) |
| MIW: atraso por ressonador (eq. 6.9) contra a cadeia simulada | — | erro < 0,1% |
| SIW: eq. 6.31 contra Deslandes e Wu | — | 0,3–0,6% |
| Linha ENZ (Fig. 6.32): Z e ε equivalente | 50 Ω, 0,3 | 48,3 Ω, 0,43 |
| Exemplo Fig. 2.34: L, C, Cls | 3,1 nH, 0,17 pF, 1,2 pF | 3,30 nH, 0,157 pF, 1,31 pF (\*\*) |

(\*) O livro obtém εeff por solver EM; aqui ele vem de fórmula fechada (εeff 6,82 contra ≈ 7,0).
(\*\*) Os valores do livro já incluem otimização e arredondamento.

## 5. Limitações conhecidas

- Microstrip: Z0 e εeff com t = 0. Perdas e dispersão são opcionais (ver seção 2) e não incluem radiação, ondas de superfície nem a perda dos patches.
- Os scripts HFSS e CST ainda modelam o cobre como condutor perfeito (PEC); σ e rugosidade não são exportados.
- Capítulo 3: a geometria (raio do CSRR, interdigital, stub) é estimativa inicial de ±20% ou mais. O livro fecha esses projetos com EM e extração (Apêndice G) ou space mapping (Apêndice H, ainda não implementado). Os modelos de circuito valem até pouco acima da ressonância dos ressonadores.
- Linhas em CPW (SRR, OSRR, OCSRR): só a extração está implementada; o projeto geométrico exige um modelo de CPW.
- Diagrama de Floquet (aba do cap. 2): o desdobramento supõe estrutura de mão direita. Linhas CRLH usam o sinal de β calculado na aba CRLH.
- Capítulo 4: dispositivos e filtros são simulados com circuitos ideais (inversores ideais, sem perdas nem parasitas). A implementação física desloca as bandas; o livro fecha com EM, extração e ASM. O α/k0 da antena leaky-wave vem de EM.
- Capítulo 6: as linhas diferenciais usam o modelo sem acoplamento magnético entre linha e CSRR (como no livro); MIW, SIW e ENZ são modelos sem perdas de condutor; a transição do SIW é estimativa inicial.
- Capítulo 5: o modelo do VLSRR usa só o acoplamento magnético (Fig. 5.8a), que o livro mostra ser dominante. A NLTL não inclui perdas do dispositivo nem da linha, e as terminações padrão são de √(L/Co); a simulação no tempo leva alguns segundos.
- A geometria da linha com carga capacitiva é **estimativa inicial**, porque não modela os efeitos de degrau e borda do patch. Em εr alto, o patch pode não caber no período: a ferramenta avisa e é preciso reduzir a "Z do patch".
- EBG com furos no terra: só a geometria e a estimativa de Bragg; a resposta exige simulação EM (eq. 2.51).
- Os scripts HFSS e CST e o `.kicad_pcb` seguem a sintaxe documentada dessas ferramentas, mas **não foram executados nas ferramentas reais** durante o desenvolvimento. Teste na versão instalada e reporte ajustes.

## 6. Material didático em Word

O arquivo `doc/Manual_LTA_Linhas_Periodicas.docx` (com PDF ao lado) é uma aula sobre o capítulo 2 do livro, com teoria, equações, capturas de tela do software, impactos de cada parâmetro e exemplos de aplicação resolvidos. Ele é gerado por script, para acompanhar a evolução do software:

```
python doc/fontes/capturar_telas.py   # abre o LTA, roda os exemplos e captura as telas (doc/fontes/img)
python doc/fontes/gerar_manual.py     # gera figuras, números e o .docx; atualiza o sumário e gera o PDF
```

Dependências extras: `python-docx`, `latex2mathml`, `Pillow`, `pymupdf` (opcional), `pywin32` e Microsoft Office instalado (conversor de equações MML2OMML e atualização do sumário). Os números dos exemplos são calculados pelo próprio núcleo do LTA (`doc/fontes/figuras.py`), então o texto acompanha mudanças nos cálculos. A captura de telas também exporta arquivos de exemplo para `saidas/exemplos_manual/`.

## 7. Estrutura do código e como estender

```
linhas_periodicas/
├── main.py, executar.bat, requirements.txt
├── lta/
│   ├── core/          cálculos (sem GUI): microstrip, abcd, carga_capacitiva, modos_acoplados,
│   │                  perfil_personalizado, floquet, aplicacoes_ebg
│   ├── geometria/     Layout (polígonos em mm) e construtores por tipo de estrutura
│   ├── exportacao/    formatos_2d (DXF, Gerber, KiCad, GDS), simuladores (HFSS, CST, Touchstone, CSV)
│   └── gui/           app.py (lista ABAS) + uma aba por arquivo
├── testes/            validação contra o livro
├── doc/               documentação
└── saidas/            arquivos exportados (padrão)
```

- **Nova ferramenta:** crie `lta/gui/aba_xxx.py` com uma classe herdando de `AbaBase` (com `TITULO` e controles em `self.painel`) e inclua-a no grupo do capítulo em `GRUPOS`, em `lta/gui/app.py`. Para um capítulo novo, acrescente um grupo. Use `estado.publicar_layout()` e `estado.publicar_resultado()` para integrar com a ferramenta Exportar.
- **Novo formato:** escreva a função em `lta/exportacao/` e registre-a em `FORMATOS` no `__init__.py`. Ela aparece automaticamente na aba Exportar.
- **Nova estrutura geométrica:** crie um construtor em `lta/geometria/construtores.py` que retorne um `Layout`. Todos os exportadores funcionam sem alteração.
