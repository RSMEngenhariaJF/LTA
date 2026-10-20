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
| **Linhas de transmissão convencionais** (cap. 0) | Microstrip, stripline, CPW, CPW com terra (GCPW), coplanar strips (CPS) e slotline: síntese Z0 → dimensões e análise; trecho de linha, stubs em aberto/curto (com desconto da ponta aberta) e transformador λ/4 (microstrip); λg, perdas, S e layout | cap. 1 (eqs. 1.68, 1.90); Cohn (1954); Simons (2001); Kirschning-Jansen-Koster (1981) |
| **Casamento de impedância** (cap. 0) | Rede L, redes Pi e T com Q especificado, stub simples (derivação e série), stub duplo, transformadores binomial e Chebyshev, tapers exponencial, triangular e Klopfenstein; limite de Bode-Fano; carta de Smith; resposta TEM ideal × microstrip e layout | Pozar, cap. 5 (Ex. 5.1–5.8); Steer 8.6 (Ex. 8.7) |
| **Linhas acopladas e acopladores** (cap. 0) | Modos par/ímpar da microstrip (Garg-Bahl) e da stripline (Cohn), síntese (Z0e, Z0o) → (W, S); acoplador λ/4, multisseção binomial e Lange; versão **LTA** com o modo ímpar carregado por capacitores periódicos (diretividade) e versão com **capacitores nas pontas**; Z0d e Z0c do par diferencial; layout de 4 portas | Pozar, 7.6–7.7 (Ex. 7.7–7.8); Martín, 2.5.2; Steer 9.10 |
| **Filtros com linhas** (cap. 0) | Passa-baixas por degraus e por stubs (Richards + Kuroda); passa-faixa de linhas acopladas (com versão **LTA** sem espúria em 2f0); rejeita-faixa e passa-faixa com stubs λ/4; ressonadores série com gaps, stubs em curto acoplados por capacitores e **combline** com capacitores | Pozar, cap. 8 (Ex. 8.5–8.10); Martín, 4.2.3.2; Steer 10.14 |
| **Ressonadores de linha** (cap. 0) | λ/2 em aberto e λ/4 em curto: comprimento (com ponta aberta), Q0 = β/2α (condutor e dielétrico), RLC equivalente; acoplamento por gap (crítico ou dado), ressonância carregada, Qe, QL e gap físico | Pozar, 6.2 e 6.6 (Ex. 6.2 e 6.6) |
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
| **Híbridas e divisores** (cap. 4) | Branch-line, rat-race, Wilkinson (igual e desigual), junção T e divisor com inversor, montados como redes de N portas; implementações CRLH dual-band, rat-race com LH, **microstrip convencional**, **slow-wave** (miniaturizada) e **concentrada** (π LC, Steer 7.13); comparação com a versão TEM ideal; faixas pelos critérios da Tabela 4.1; layouts | 4.2.1–4.2.2 (eqs. 4.2–4.19, Figs. 4.5 e 4.16); 2.5.2; Pozar 7.2–7.5 e 7.8 |
| **Filtros passa-faixa e dual-band** (cap. 4) | Protótipos Butterworth/Chebyshev; passa-faixa com inversores e ressonadores; transformação dual-band; filtro periódico com células híbridas de CSRR | 4.2.3 (eqs. 4.38–4.50) |
| **Linhas sintonizáveis** (cap. 5) | Varactor, HBV e RF-MEMS (pull-in); linha com SRR e varactor (rejeita-faixa sintonizável); CRLH balanceada sintonizável como defasador (FoM) e antena leaky-wave de varredura em frequência fixa | 5.2–5.3 (eqs. 5.1–5.12) |
| **Linha não linear (NLTL)** (cap. 5) | Sólitons (atraso × amplitude, forma), decomposição de pulsos e multiplicação de frequência por simulação no tempo | 5.4 (eqs. 5.13–5.25) |
| **Linha diferencial (filtro de modo comum)** (cap. 6) | Par diferencial com CSRR/DS-CSRR: zero do modo comum, banda máxima de rejeição (estrutura infinita), filtro de N células para os modos comum e diferencial, banda × ordem | 6.3 (Fig. 6.12, eq. 6.11, Tabela 6.1; eqs. 3.93–3.97) |
| **Linhas de atraso MIW/EIW** (cap. 6) | Cadeias de SRR/CSRR acoplados: dispersão, banda, fator de onda lenta; etiqueta RFID por reflectometria (bits, ressonadores) e ecos de um pulso | 6.2 (eqs. 6.1–6.10) |
| **Guia integrado ao substrato (SIW)** (cap. 6) | Largura, vias e banda monomodo pelo corte do TE10; dispersão, λg, Z_TE, perda dielétrica; layout com vias e transições em cunha | 6.5 (eqs. 6.31–6.34) |
| **Linha com capacitores negativos** (cap. 6) | Linha ENZ não Foster: C' efetiva, estabilidade, Z, ε equivalente, vp e vg superluminais | 6.4.2 (Fig. 6.32) |
| **Conectores SMA e transições** (prático) | Transição coaxial-microstrip como Ls série + Cs derivação (predefinida, informada ou pelo VSWR); trecho de compensação junto ao pad; fixture coaxial-linha-coaxial; de-embedding de .s2p | Garg-Bahl-Bianchi 1.4.1 e Tabela 2.5 |
| **Caracterização por medidas** (prático) | Duas linhas (sem os conectores), uma linha (Zc, γ, RLGC), dois ressonadores, Q0 → α, VSWR → Z0; εr pela inversão de Hammerstad-Jensen; calibração **TRL** e renormalização de S; exemplos sintéticos | Garg-Bahl-Bianchi 1.5 e 2.4.7; Steer 6.10-6.11 |
| **Tolerâncias, potência e orçamento de erros** (prático) | Sensibilidades, pior caso e Monte Carlo (microstrip e acopladas), ΔW exigido; potência média (térmica) e de pico (ruptura); orçamento de VSWR como na Tabela 2.5 | Garg-Bahl-Bianchi 2.4.2-2.4.6 e 4.5.3 |
| **DGS** (prático) | Haltere como LC paralelo (eq. 6.8), linha DGS periódica (rejeição e onda lenta), passa-baixas com DGS e stubs, extração do modelo π de .s2p; layout com as fendas no terra | Garg-Bahl-Bianchi cap. 6 |
| **Balun de Marchand** (prático) | Duas seções λ/4 acopladas, acoplamento ótimo (C = 1/√(1 + ZL/Zs)), síntese em microstrip, versão **LTA** com fase de 180° e versão **miniaturizada** (seções θ0 < 90° + capacitores ao terra); layout com vias | Garg-Bahl-Bianchi 4.9; Steer 13.2 |
| **Amplificador: estabilidade, ganho e casamento** (Steer) | S do pHEMT do livro ou de .s2p; k, \|Δ\|, μ; círculos de estabilidade na Smith; GTU,max, GMA, GMS e U; casamento conjugado simultâneo com rede L ou stub e GT resultante | Steer 11.5-11.7 (Fig. 11-8, Tabelas 11-6 e 11-7) |
| **Atenuadores resistivos** (Steer) | Pads T e Pi com Z01 ≠ Z02, atenuação mínima, versões balanceadas, valores E24/E96 e potência por resistor | Steer 7.3.2 (eqs. 7.2-7.14) |
| **Exportar** | Layout para KiCad, DXF, Gerber, GDSII, script HFSS e macro CST; parâmetros S em Touchstone e CSV | — |

Cada aba de projeto tem botões **Exemplo** que reproduzem figuras do livro, **Enviar p/ exportação** (manda o layout para a aba Exportar) e **Exportar S2P**.

### Capítulo 0: linhas de transmissão convencionais

A primeira ferramenta do menu projeta as linhas planares comuns, que são a base das linhas artificiais ([linhas_planares.py](../lta/core/linhas_planares.py)).

- **Modelos (quase-estáticos, metal de espessura nula):**
  - microstrip de Hammerstad e Jensen, com dispersão e perdas;
  - stripline centrada de Cohn, exata por transformação conforme;
  - CPW, GCPW e CPS em substrato de espessura h, de Simons (integrais elípticas).
- **Síntese:** a largura (W, ou a fita central S, ou as fitas da CPS) é obtida por bisseção, com a fenda G (ou a separação S da CPS) e a distância b da stripline fixas.
- **Elementos:**
  - trecho de linha com o comprimento elétrico pedido em f0;
  - stub em aberto ou em curto em derivação: Yin = tanh(γl)/Z0 ou 1/(Z0·tanh γl);
  - transformador λ/4 com Zt = √(Z0·Z2), simulado com portas de impedâncias diferentes.
  - Stubs e transformador estão disponíveis em microstrip.
- **Perdas:** dielétrico pela eq. 1.68 em todas as linhas quase-TEM; condutor (eq. 1.90/Pucel) só na microstrip.
- **Layout e exportação:**
  - microstrip (trecho, stub com pad e via, transformador em degraus);
  - CPW e CPS sem plano de terra inferior: os arquivos 2D omitem a camada TERRA, e no HFSS/CST o terra fica inerte. As portas de onda precisam de ajuste, como indicado no LEIA-ME;
  - GCPW com cerca de vias opcional.
  - A stripline (três camadas) tem só análise.
- **Pré-visualização:** o corte transversal passa a mostrar todo o cobre do topo que cruza a seção (planos coplanares, fitas) e as vias.
- **Ponta aberta:** o stub em aberto tem o comprimento físico descontado de Δl (Kirschning-Jansen-Koster); a resposta usa sempre o comprimento físico + Δl, de modo que desmarcar o desconto mostra o deslocamento da ressonância.

### Capítulo 0: circuitos com linhas (Pozar)

Quatro abas aplicam o livro de D. M. Pozar (*Microwave Engineering*, 4ª ed., 2012) às linhas convencionais. Todas comparam a rede TEM ideal (tracejado) com a microstrip no substrato atual (contínuo), com perdas opcionais, mostram os limites do substrato e geram layout.

- **Correções físicas** ([descontinuidades.py](../lta/core/descontinuidades.py)):
  - extremidade aberta Δl (Kirschning, Jansen e Koster, 1981), usada em stubs, ressonadores e filtros;
  - gap série em modelo π (Garg e Bahl, 1978) e síntese da largura do gap para uma capacitância;
  - limites do substrato: TE1 e TM1 de onda de superfície (Pozar, eq. 3.174) e modo transversal da microstrip, com aviso quando f máx se aproxima deles;
  - fórmulas fechadas do Pozar (eqs. 3.195–3.199), para comparação na aba Substrato.
- **Casamento** ([casamento.py](../lta/core/casamento.py)):
  - carga R + jX em f0, modelada como R em série com L ou C (como nos exemplos do livro);
  - rede L com as duas topologias (eqs. 5.3 e 5.6), stubs simples e duplo por busca das raízes de Re{Y} = Y0 (ou Re{Z} = Z0), comprimentos em λ e em mm;
  - binomial (eq. 5.53) e Chebyshev pela expansão de T_N(sec θm·cos θ) em cossenos (eqs. 5.61–5.64), válido para qualquer N;
  - tapers com o perfil Z(z) do livro (Klopfenstein pela integral φ(x, A) com I1 em série) e resposta exata por cascata de trechos curtos;
  - Bode-Fano para RC e RL série (Fig. 5.23).
- **Linhas acopladas** ([linhas_acopladas.py](../lta/core/linhas_acopladas.py)):
  - a rede de 4 portas sai dos ABCD dos modos par e ímpar; cada modo tem Z, εeff e perdas próprios, então a diretividade finita da microstrip aparece naturalmente;
  - multisseção binomial com N ímpar qualquer: as derivadas pares de C(θ) em π/2 são anuladas por sistema linear (Ex. 7.8 para N = 3);
  - Lange: Z0e e Z0o do par pela eq. 7.100 e modos de 4 linhas pela eq. 7.97 (layout não gerado por causa dos fios de solda).
- **Acopladores com LTA:** na microstrip, εeff do modo par é maior que o do ímpar. Capacitores Cc entre as linhas, em células periódicas, só carregam o modo ímpar (2Cc para o plano de simetria). A geometria é sintetizada para Zo' = Zo·√(εee/εeo) e Cc baixa Zo' ao alvo enquanto iguala as velocidades. Isso é a célula slow-wave do Martín (seção 2.5.2) aplicada a um modo. Em RO4003C (Kirschning-Jansen), um acoplador de 10 dB passa de 13 dB para 70 dB de diretividade; o de 3 seções e 20 dB, de 1 para 48 dB. Seções fracas demais para compensar ficam convencionais (aviso no texto).
- **Filtros** ([filtros_convencionais.py](../lta/core/filtros_convencionais.py)):
  - degraus de impedância (eq. 8.86) com comparação ao protótipo concentrado;
  - Richards + Kuroda com linhas redundantes nas pontas (N ≤ 3 começando com L, ou N ≤ 2 com C; acima disso, use degraus);
  - linhas acopladas: Z0·J (eq. 8.121), Z0e e Z0o (eq. 8.108), W e S sintetizados e resposta com as duas pontas abertas de cada seção; a versão **LTA** suprime a espúria em 2f0 (cerca de −80 dB contra 0 dB na microstrip comum), como discutido em Martín 4.2.3.2;
  - stubs λ/4 (eqs. 8.130–8.131), ressonadores com gaps (eqs. 8.134–8.135, gap físico por Garg-Bahl e Cp descontado do ressonador) e stubs em curto com capacitores (eqs. 8.136–8.141).
- **Ressonadores** ([ressonadores.py](../lta/core/ressonadores.py)): Q0 com as perdas do modelo da aba Substrato e comparação com αc ≈ Rs/(Z0W) do Pozar; acoplamento crítico bc = √(π/(2Q0)) (λ/2) ou √(π/(4Q0)) (λ/4); ressonância deslocada pela raiz de tan βℓ + bc = 0; Q de radiação (Garg-Bahl eq. 2.142) informado à parte.

### Projeto prático (Garg, Bahl e Bianchi)

Do livro *Microstrip Lines and Slotlines* (4ª ed., Artech House, 2024) vieram modelos mais precisos para as linhas existentes e um grupo novo de ferramentas.

- **Microstrip com espessura** ([microstrip.py](../lta/core/microstrip.py)): opção do substrato ("Espessura t em Z0 e εeff"); largura efetiva We (eq. 2.125) na impedância em ar de Hammerstad-Jensen e εeff − (εr − 1)/4,6·(t/h)/√(W/h) (eqs. 2.126-2.127). Desligada por padrão, para manter as validações do Martín com t = 0. Afeta todas as abas.
- **Coplanares e slotline** ([linhas_planares.py](../lta/core/linhas_planares.py)): espessura em CPW, GCPW e CPS pelas larguras efetivas do semiplano superior (eqs. 7.110-7.114); perda no condutor das coplanares (eq. 7.126); slotline dispersiva de Janaswamy-Schaubert (2,22 ≤ εr ≤ 9,8) e Garg-Gupta (εr até 20), com síntese em f0 e layout sem terra.
- **Descontinuidades** ([descontinuidades.py](../lta/core/descontinuidades.py)): degrau de largura (Δl = Δloc·(1 − W2/W1)), curva de 90° (Cb e Lb, eqs. 3.57-3.58), chanfro ótimo M = 52 + 65·e^(−1,35W/h) (eq. 3.115) com o encurtamento Δb (eq. 3.116), via de aterramento (L e R de Goldfarb-Pucel, eqs. 3.61-3.62), entalhe (eq. 3.59), fT de acoplamento com a onda de superfície (eq. 2.140) e Q de radiação (eq. 2.142).
  - Os stubs em curto agora terminam na via real em todas as abas, e o comprimento é descontado do equivalente da via quando a correção está marcada.
  - O filtro de degraus e os transformadores multisseção corrigem os comprimentos pelos degraus, inclusive com as alimentações, e o modelo microstrip inclui esses degraus.
  - Os layouts do Wilkinson, da junção T, do rat-race, do acoplador e do balun podem ter as curvas chanfradas; no Wilkinson, cada braço é alongado de 2·Δb.
- **Linhas acopladas: Kirschning-Jansen** ([linhas_acopladas.py](../lta/core/linhas_acopladas.py)), agora o modelo padrão:
  - εeff estático (eqs. 4.96-4.98, erro ≤ 0,7%) e impedâncias (eqs. 4.100-4.101, ≤ 0,6%); Garg-Bahl continua selecionável;
  - dispersão de cada modo (eqs. 4.79, 4.81-4.82), opcional;
  - espessura pelas capacitâncias: ΔCf da linha isolada e Cgt = 2ε0t/S no modo ímpar (eqs. 4.108-4.109);
  - perda no condutor de cada modo pela indutância incremental de Wheeler (eqs. 4.110-4.113), calculada numericamente recuando todas as paredes metálicas; o modo ímpar perde mais com linhas próximas;
  - ponta aberta de cada modo (eq. 4.155), sempre incluída no filtro de linhas acopladas como extensão elétrica Δle/Δlo.
- **LTA com ponta aberta:** no filtro de linhas acopladas, Δle ≠ Δlo volta a desigualar as fases dos modos. A compensação LTA passou a igualar o comprimento elétrico total, (ℓ + Δlo)·√εeo' = (ℓ + Δle)·√εee. Com isso a espúria em 2f0 cai para cerca de −80 dB.
- **Conectores** ([conectores.py](../lta/core/conectores.py)):
  - transição coaxial-microstrip como Ls série + Cs derivação (Fig. 1.20);
  - predefinidos que reproduzem o VSWR da Tabela 2.5 (1,065 em 3 GHz e 1,10 em 10 GHz, Cs ≈ 30 fF): são ordens de grandeza, não dados de fabricante;
  - compensação por um trecho curto de alta impedância (ou baixa, se a transição for indutiva) com θ/ω·(Zc − Z0²/Zc) = Z0²·(Cs − Ls/Z0²);
  - de-embedding: M_DUT = A⁻¹·M·B⁻¹, com A e B transição + trecho + alimentação.
- **Caracterização** ([caracterizacao.py](../lta/core/caracterizacao.py)):
  - duas linhas: os autovalores de M1⁻¹·M2 dão e^(±γΔL) sem a influência dos conectores;
  - uma linha: Zc e γ pelas eqs. 2.147-2.148 e RLGC;
  - dois ressonadores (eqs. 1.94-1.95), Q0 → α (eq. 1.96) e VSWR → Z0 (eq. 1.91);
  - εr pela inversão de Hammerstad-Jensen (mais precisa que a fórmula F do livro);
  - a fase desenrolada é ancorada para que β ∝ f passe pela origem.
- **Robustez** ([robustez.py](../lta/core/robustez.py)):
  - sensibilidades numéricas S = (B/A)·∂A/∂B sobre os mesmos modelos do programa; pior caso (eqs. 2.97-2.98) e ΔW exigido (eqs. 2.100-2.101);
  - Monte Carlo uniforme; acopladas (eq. 4.125);
  - potência média pelo modelo térmico de placas paralelas (eqs. 2.84-2.95), com K da Tabela 2.2 ou de datasheet (laminados Rogers: conferir com o fornecedor);
  - potência de pico pela ruptura do dielétrico (eq. 2.96);
  - orçamento de VSWR como a Tabela 2.5, com o pior caso combinado pelo produto dos VSWR.
- **DGS** ([dgs.py](../lta/core/dgs.py)):
  - C = ωc/(2Z0(ω0² − ωc²)) e L = 1/(ω0²C) a partir da rejeição e do corte de 3 dB da simulação EM (eq. 6.8);
  - R pelas perdas (eq. 6.9) e extração do modelo π (eqs. 6.10-6.13) de um .s2p;
  - linha DGS periódica analisada como célula ABCD (Bloch), dando a banda de rejeição e o fator de onda lenta;
  - passa-baixas de 3 polos: DGS como indutores série e um par de stubs abertos como capacitor; o exemplo de 1,3 GHz dá 3,85 nH e 2,38 pF, contra 3,678 nH e 2,33 pF do livro, que vêm de otimização;
  - as dimensões do haltere só desenham o layout e devem ser ajustadas em EM.
- **Balun de Marchand** ([balun.py](../lta/core/balun.py)):
  - duas seções acopladas montadas na rede nodal (novo `Rede.multiporta`), linhas b aterradas nas pontas externas;
  - o acoplamento ótimo em TEM reproduz C = 1/√(1 + ZL/Zs);
  - acoplamentos fortes (ZL ≤ 100 Ω com Zs = 50 Ω) não cabem em linhas lado a lado, e a ferramenta avisa;
  - em microstrip, a fase entre as saídas cai para −173°; a versão LTA (modo ímpar compensado) volta a 180,0°.

### Sistemas RF (Steer)

Do livro de M. Steer, *Microwave and RF Design: A Systems Approach* (SciTech, 2010), vieram um grupo novo e opções nas abas existentes.

- **Calibração TRL e renormalização** ([caracterizacao.py](../lta/core/caracterizacao.py), aba Caracterização):
  - TRL de Engen-Hoer com thru de comprimento zero: os autovetores de Line·Thru⁻¹ (matrizes T) dão as caixas de erro a menos de escalas, e o Reflect (desconhecido, igual nas duas portas) fixa a razão; só o sinal do Reflect (aberto ou curto) precisa ser informado;
  - o exemplo sintético gera Thru, Line (λ/4 no centro da banda), Reflect e DUT com transições SMA e confere o DUT corrigido com o verdadeiro (erro ~10⁻¹⁵);
  - avisa quando a fase da Line sai de 20°-160° (TRL mal condicionado);
  - renormalização de S para outra Zref via matriz Z (Steer 6.10), por porta.
- **Amplificador** ([amplificador.py](../lta/core/amplificador.py), aba nova):
  - k de Rollett, |Δ|, μ e μ' de Edwards-Sinsky; círculos de estabilidade de fonte e carga;
  - GTU,max, GMA (k > 1), GMS = |S21/S12| e U de Mason;
  - casamento conjugado simultâneo e as redes de entrada e saída (rede L ou stub do Pozar), com o GT na banda;
  - S do pHEMT FPD6836P70 da Fig. 11-8 embutido; .s2p em outra Zref é renormalizado.
- **Redes Pi e T com Q** ([casamento.py](../lta/core/casamento.py), aba Casamento): duas redes L através de RV = Rmaior/(Q² + 1) (Pi) ou Rmenor·(Q² + 1) (T); as quatro soluções, como na Fig. 8-30, com o elemento do meio somado.
- **Híbridas concentradas** ([dispositivos.py](../lta/core/dispositivos.py), aba Híbridas e divisores, implementação "Concentrada"): cada linha vira o π LC equivalente em f0 (série X = Z·sen θ, derivações B = (1 − cos θ)/(Z·sen θ)); 90° é passa-baixas (L = Z/ω0, C = 1/(ω0Z)) e 270° vira −90° passa-altas. Sem layout.
- **Acoplador com capacitores nas pontas** ([linhas_acopladas.py](../lta/core/linhas_acopladas.py), aba Linhas acopladas): o modo ímpar, mais rápido, tem θx = 90°·√(εeo/εee); a linha é ressintetizada com Zx = Zo/sen θx e um capacitor c = cos θx/(2ω0·Zx·sen θx) entre as linhas em cada ponta (Steer 9.10). A diretividade em f0 passa de 60 dB, mas cai para ~30 dB a ±25%; a LTA distribuída mantém > 60 dB numa oitava.
- **Impedâncias de modo diferencial e comum** (aba Linhas acopladas, análise): Z0d = 2·Z0o e Z0c = Z0e/2.
- **Combline com capacitores** ([combline.py](../lta/core/combline.py), aba Filtros com linhas): procedimento da seção 10.14 generalizado para ordem N (protótipo → inversores → ressonadores LC → capacitor + stub em curto de 45° em f0 → escala dos nós → inversores em π de stubs → nível interno → inversores de entrada por Cb série); análise por rede nodal com stubs TEM e layout com vias embaixo, capacitores no topo e Cb nas alimentações (W e S estimados por pares, a ajustar em EM).
- **Marchand miniaturizado** ([balun.py](../lta/core/balun.py), aba Balun): seções de θ0 < 90° com capacitores ao terra no meio da linha a, na ponta aberta e nas saídas; Nelder-Mead ([otimizacao.py](../lta/core/otimizacao.py)) primeiro com o acoplamento convencional e depois refinando também o acoplamento (fica 0,2-1 dB mais forte); depois resintetiza em microstrip e reajusta os capacitores.
- **Atenuadores** ([passivos.py](../lta/core/passivos.py), aba nova): pads T e Pi (eqs. 7.2-7.11), atenuação mínima com Z01 ≠ Z02 (eq. 7.14), versões balanceadas, valores E24/E96 com o erro resultante e potência por resistor.

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
  - **Rede de N portas:** cada linha entra pela sua matriz ABCD na matriz Y nodal. As portas são terminadas nas impedâncias de referência (podem ser diferentes por porta) e o sistema completo é invertido. Antes, os nós internos eram eliminados por redução de Kron, que ficava singular quando um nó só via linhas λ/4 ideais exatamente em f0.
  - **Dispositivos:**
    - branch-line: entrada 1, direta 2, acoplada 3, isolada 4;
    - rat-race: porta 4 = Σ, porta 1 = Δ, como na Fig. 4.7;
    - Wilkinson igual e desigual (P3/P2 = K², eq. 7.37 do Pozar), com ou sem transformadores λ/4 de saída;
    - junção T sem perdas (Pozar, Ex. 7.1), com portas de referência Z2/Z3 ou transformadores;
    - divisor em Y com inversor (Figs. 4.11 e 4.15).
  - **Linhas:**
    - CRLH dual-band (−90° em f1, +90° em f2);
    - PLH de −90° no lugar do ramo de 270° do rat-race (banda larga, Figs. 4.1 e 4.5);
    - microstrip convencional (Pozar), com dispersão e perdas opcionais;
    - **slow-wave:** linha não carregada de Zo com carga capacitiva periódica (stubs abertos ou patches). A célula distribuída é resolvida exatamente em f0 para dar a fase φ/N e a impedância de Bloch Z pedidas; o fator de redução é ≈ Z/Zo. Na branch-line de 2 GHz em RO4003C, com Zo = 100 Ω e 3 células, os ramos caem para 39% e 53% do comprimento, com banda (|S11| e isolamento < −15 dB) de 16% contra 18% da convencional. Stubs voltados para fora do quadrado evitam a sobreposição que patches largos causam nos cantos.
  - **Comparação:** a versão escolhida contra a mesma rede com linhas TEM ideais em f0.
  - **Layouts:** branch-line (convencional e slow-wave), rat-race em anel, Wilkinson em "C" com gap para resistor SMD e junção T. As CRLH e LH são concentradas e ficam sem layout.
  - **Métricas:** faixa contínua em torno de cada frequência de projeto em que valem |S11| < −10 dB, saídas no nível esperado ± 1 dB (−3 dB, ou a divisão K² nos desiguais), balanço de fase ± 5° e isolamento < −15 dB.
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
| Stripline 50 Ω, b = 3,16 mm, εr = 2,2 (Pozar, ex. 3.5) | W ≈ 2,66 mm (fórmula aproximada) | 2,62 mm (Cohn, exata) |
| CPW no ar, k = 0,5 | 30π·K(k')/K(k) = 120,6 Ω | 120,6 Ω |
| CPW em substrato muito espesso: εeff | (εr + 1)/2 | idem (10⁻⁵) |
| GCPW com h → ∞ | igual à CPW | igual (10⁻⁴) |
| CPW × CPS duais (Booker): Z1·Z2 | η0²/(4εeff) | a 0,14% |
| Stub λ/4 aberto / em curto; transformador λ/4 em f0 | zero / passagem / casado | < −60 dB / 0 dB / < −60 dB |
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

**Pozar, *Microwave Engineering*, 4ª ed. (classes `Pozar*` e `LayoutsPozar` dos testes):**

| Item | Livro | Ferramenta |
|---|---|---|
| Ex. 3.7: microstrip 50 Ω em alumina 0,5 mm | W = 0,483 mm, εe = 6,665 | 0,4828 mm, 6,664 |
| Ex. 5.1: rede L, ZL = 200 − j100 Ω, 500 MHz | b = 0,29, x = 1,22; 0,92 pF / 38,8 nH; 2,61 pF / 46,1 nH | 0,290, 1,225; 0,923 pF / 39,0 nH; 2,60 pF / 46,1 nH |
| Ex. 5.2: stub em curto em derivação | d = 0,110λ / 0,260λ; ℓ = 0,095λ / 0,405λ | 0,110 / 0,259; 0,095 / 0,405 |
| Ex. 5.3: stub aberto em série | d = 0,120λ / 0,463λ; ℓ = 0,397λ / 0,103λ | 0,120 / 0,463; 0,398 / 0,102 |
| Ex. 5.4: stub duplo aberto, λ/8 | b1 = 1,314 / −0,114; b2 = 3,38 / −1,38; ℓ = 0,146, 0,204 / 0,482, 0,350λ | idem |
| Ex. 5.6: binomial, 3 seções, 100 → 50 Ω | 91,7 / 70,7 / 54,5 Ω; banda 70% | 91,70 / 70,71 / 54,53 Ω; 70,3% |
| Ex. 5.7 e Tabela 5.2: Chebyshev, Γm = 0,05 | 57,5 / 70,7 / 87,0 Ω; θm = 44,7° | 57,48 / 70,71 / 86,99 Ω; 44,73° (N = 2 e 4: a 0,3% da tabela exata) |
| Ex. 5.8: Klopfenstein, Γm = 0,02 | A = 3,543; \|Γ\| ≤ Γm para βL ≥ A | 3,545; \|Γ\| ≤ 0,024 (cascata exata) |
| Ex. 6.2: λ/2 em Teflon, 5 GHz | W = 5,08 mm, εe = 1,80, ℓ = 2,24 cm, αd = 0,024 Np/m, αc = 0,0724 Np/m | 5,078 mm, 1,797, 2,236 cm, 0,0240, 0,0723 (fórmula do Pozar) |
| Ex. 6.2: Q0 | 783 (com β = 151 rad/m) | 729: β = 2πf√εe/c = 140,5 rad/m com os dados do próprio exemplo (\*\*\*) |
| Ex. 6.6: gap crítico | bc = 0,05, C = 0,032 pF, f = 4,918 GHz | 0,0500, 0,0319 pF, 4,920 GHz |
| Ex. 7.1: junção T 2:1 | Z = 150 e 75 Ω; Γ = −0,666 e −0,333 | idem |
| Eq. 7.37: Wilkinson desigual P3/P2 = 2 | divisão 1/3–2/3, casado e isolado | exato (com e sem transformadores) |
| Ex. 7.7: acoplador 20 dB em stripline, b = 3,2 mm | Z0e/Z0o = 55,28/45,23 Ω; W = 2,59 mm, S = 0,98 mm | 55,28/45,23 Ω; 2,597 mm, 0,980 mm |
| Ex. 7.8: binomial de 3 seções | C = 0,0125/0,125; 50,63/49,38 e 56,69/44,10 Ω | idem; 19,97 dB em f0 (teoria de acoplamento fraco: 20 dB) |
| Fig. 7.30: microstrip acoplada, εr = 10, W/h = 1, S/h = 0,5 | ≈ 60 / 38 Ω (gráfico) | 59,9 / 37,7 Ω (Garg-Bahl) |
| Ex. 8.5: Richards + Kuroda, 3 dB, N = 3 | 1,299 / 4,350 / 1,405; 217,5, 64,9, 70,3 Ω | idem; identidade exata (10⁻¹⁵) |
| Ex. 8.6: degraus, N = 6 | βℓ = 11,8 / 33,8 / 44,3 / 46,1 / 32,4 / 12,3°; W = 11,3 / 0,428 mm | 11,9 / 33,8 / 44,3 / 46,1 / 32,4 / 12,4°; 11,3 / 0,428 mm |
| Ex. 8.7: linhas acopladas, 0,5 dB, N = 3 | Z0J = 0,3137 / 0,1187; Z0e/Z0o = 70,61/39,24, 56,64/44,77; ~20 dB em 1,8 GHz | idem; 19,4 dB |
| Ex. 8.8: stubs λ/4, Δ = 15% | 265,9 / 387,0 Ω | idem |
| Ex. 8.9: gaps | C = 0,554 / 0,192 pF; θ = 155,8 / 166,5° | idem |
| Ex. 8.10: stubs em curto + capacitores | C = 0,2896 / 0,0756 pF; ℓ = 73,6 / 83,2° | idem |
| LTA: acoplador 10 dB em RO4003C | — | diretividade 13 → 70 dB |
| LTA: filtro de linhas acopladas, espúria em 2f0 | — | 0 dB → cerca de −80 dB |

**Garg, Bahl e Bianchi, *Microstrip Lines and Slotlines*, 4ª ed. (classes `GargBahl*` dos testes):**

| Item | Livro | Ferramenta |
|---|---|---|
| Fig. 4.12: acopladas, εr = 9,7, W/h = 1; εeff estático | par ≈ 7,0; ímpar 5,55 (S/h = 0,1) a 5,8 (S/h = 1) | 7,00; 5,53 a 5,82 |
| Fig. 4.12: εeff em 30 GHz | par ≈ 8,35; ímpar ≈ 6,05 (S/h = 0,1) e 6,85 (S/h = 1) | 8,33; 6,03 e 6,85 |
| Garg-Bahl × Kirschning-Jansen | ~3% (Garg-Bahl) | 1-2% em S/h = 0,5; até 6% em Zo com S/h = 0,1 |
| Wheeler: linhas afastadas | perda = da linha isolada | a 5%; αco > αce com linhas próximas |
| Ponta aberta acoplada, S grande | Δle → Δloc | idem (10⁻³ mm) |
| Via h = 380 µm, r = 100 µm | 67 pH | 67,5 pH |
| fT, quartzo h = 0,5 mm | > 100 GHz | 106 GHz |
| Qr, εr = 3,8, 50 Ω | 0,4252·10⁴/(fh)² | 4246 em 1 GHz·mm |
| Chanfro ótimo, alumina 50 Ω | ~72% | 69% (eq. 3.115) |
| Espessura, alumina t/h = 0,016 | Δεeff −0,03 a −0,08 (Tabela 2.5) | −0,030 |
| Perda no condutor da CPW (Fig. 7.23) | 0,45-0,75 dB/cm (5-20 GHz) | 0,54-0,85 dB/cm (eq. 7.126 fora de t > 3δ) |
| Potência média, alumina 50 Ω em 2/10/20 GHz | 12,1 / 5,2 / 3,4 kW (Tabela 2.3) | 14,2 / 6,2 / 4,2 kW |
| DGS f0 = 8, fc = 3,5 GHz | — | L = 3,677 nH, C = 0,1076 pF; −3,01 dB em fc |
| Passa-baixas DGS 1,3 GHz | Ld = 3,678 nH, Cs = 2,33 pF (otimizados) | 3,85 nH, 2,38 pF (protótipo 0,01 dB) |
| Marchand, C ótimo | — | 1/√(1 + ZL/Zs) a 10⁻⁴ |
| Transição com VSWR 1,10 em 10 GHz | Tabela 2.5 | Cs = 30,3 fF; compensação: perda de retorno 27,6 → 37 dB |
| Duas linhas e de-embedding (dados sintéticos com transições) | — | εeff, α e S exatos (10⁻⁹) |
| TRL sintético (aberto e curto) e renormalização | Steer 6.10-6.11 | DUT e ΓR exatos (10⁻⁹); ida e volta 10⁻¹² |
| pHEMT a 0,5 GHz: k, \|Δ\| | 0,15178 / 0,62757 (Tabela 11-6) | 0,15178 / 0,62757 |
| GTU,max 0,5-5 GHz | 36,62 / 31,07 / 24,88 / 21,26 / 18,73 / 16,00 dB (Tabela 11-7) | idem, mas 17,00 dB em 5 GHz (\*\*\*\*) |
| Círculo de fonte a 0,5 GHz | centro 0,92820 + j0,50682, raio 0,22434 | idem (10⁻⁴) |
| Casamento conjugado a 8 GHz | GT = GMA | 13,95 dB |
| Combline 1 GHz, 10%, ε = 0,1 | C1 = 1,6515, C2 = 2,7171 pF; Z1 = 67,75 Ω; 442,38 Ω; Cb = 1,22306 pF | 1,6428 / 2,7076 pF; 67,77 Ω; 443,24 Ω; 1,2217 pF (0,6%) |
| Rede Pi, Ex. 8.7 (50 → 500 Ω, Q = 10) | X1 = 16,4; Xa = 14,933; Xb = 49,5; X3 = 50 Ω | X1 = 16,58 (\*\*\*\*\*); 14,934; 49,505; 50 Ω |
| Pads de 10 dB em 50 Ω | T: 25,97 / 35,14 Ω; Pi: 96,25 / 71,15 Ω | idem |
| Pad de perda mínima 50 → 75 Ω (eq. 7.14) | 5,72 dB | 5,719 dB |
| Acoplador 10 dB em RO3010 com capacitores nas pontas | — | diretividade em f0: 12 → 131 dB |
| Marchand miniaturizado, θ0 = 45°, 50:200 Ω | — | \|S11\| < −100 dB, fase 180,0° em f0 |

(\*\*\*\*) Com as próprias S da Fig. 11-8, GTU,max em 5 GHz é 17,00 dB; o 16,00 da Tabela 11-7 parece erro de digitação (os vizinhos conferem).

(\*\*\*\*\*) X1 = Rs/Qesq com Qesq = √(50/4,9505 − 1) = 3,0166 dá 16,58 Ω; o livro imprime 16,4 Ω, mas os demais valores (14,933 e 64,433) são coerentes com 16,58.

(\*\*\*) O livro escreve β = 151,0 rad/m, mas 2π·5 GHz·√1,80/c = 140,5 rad/m; com os αc e αd do próprio exemplo, Q0 = 140,5/(2·0,0964) ≈ 729.

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
- Capítulo 0 (Pozar): junções T, curvas, cantos e degraus não são compensados (apenas a ponta aberta e, no filtro com gaps, a capacitância Cp). As perdas dos modos das linhas acopladas usam a linha isolada como aproximação. O gap de Garg-Bahl vale para 0,5 ≤ W/h ≤ 2 e 0,1 ≤ s/W ≤ 1; fora disso é extrapolado. Os capacitores de compensação dos acopladores e filtros LTA não são desenhados: as posições ficam como cotas e no LEIA-ME, para SMD 0201 ou dedos interdigitais a ajustar em EM. O Lange não tem layout por causa dos fios de solda.
- Garg-Bahl: a junção T e o cruzamento só têm dados em gráficos no livro e não foram modelados; a transição SMA usa Ls e Cs informados ou estimados pelo VSWR, não dimensões de conector; a slotline cobre W/λ0 ≤ 0,075 (fendas estreitas) e tem só perda dielétrica; a potência não inclui resistências térmicas de solda e fixação; as dimensões do haltere da DGS são só desenho; o balun desconsidera vias, junções e a separação entre as seções.
- Steer: o combline é analisado com stubs TEM ideais; W e S da microstrip vêm da aproximação por pares e o meio não é homogêneo, então o layout é ponto de partida para EM. Os capacitores do combline, do Marchand miniaturizado e do acoplador com capacitores nas pontas são marcados como cotas (SMD), sem parasitas. O amplificador não inclui polarização, ruído nem estabilização; as redes de casamento são ideais. O TRL supõe thru de comprimento zero e Line com a mesma Z0 do Thru.
- Slow-wave: o patch modela só a capacitância excedente (sem a indutância e o degrau); o stub aberto usa a admitância real tan(βℓ)/Z, mas sem a junção T. Os layouts slow-wave são gerados só para a branch-line.
- Os scripts HFSS e CST e o `.kicad_pcb` seguem a sintaxe documentada dessas ferramentas, mas **não foram executados nas ferramentas reais** durante o desenvolvimento. Teste na versão instalada e reporte ajustes.

## 6. Material didático em Word

O arquivo `doc/Manual_LTA_Linhas_Periodicas.docx` (com PDF ao lado) é uma aula sobre os capítulos 1 a 6 do livro (linhas convencionais, periódicas, metamateriais, aplicações, sintonizáveis e não lineares e outras linhas avançadas), com teoria, equações, capturas de tela do software, impactos de cada parâmetro e exemplos de aplicação resolvidos. Ele é gerado por script, para acompanhar a evolução do software:

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
│   │                  perfil_personalizado, floquet, aplicacoes_ebg, ... e, do Pozar:
│   │                  descontinuidades, casamento, linhas_acopladas, filtros_convencionais, ressonadores
│   │                  e, de Garg-Bahl-Bianchi: conectores, caracterizacao, robustez, dgs, balun
│   ├── geometria/     Layout (polígonos em mm), construtores (linhas artificiais) e circuitos (Pozar)
│   ├── exportacao/    formatos_2d (DXF, Gerber, KiCad, GDS), simuladores (HFSS, CST, Touchstone, CSV)
│   └── gui/           app.py (lista ABAS) + uma aba por arquivo
├── testes/            validação contra o livro
├── doc/               documentação
└── saidas/            arquivos exportados (padrão)
```

- **Nova ferramenta:** crie `lta/gui/aba_xxx.py` com uma classe herdando de `AbaBase` (com `TITULO` e controles em `self.painel`) e inclua-a no grupo do capítulo em `GRUPOS`, em `lta/gui/app.py`. Para um capítulo novo, acrescente um grupo. Use `estado.publicar_layout()` e `estado.publicar_resultado()` para integrar com a ferramenta Exportar.
- **Novo formato:** escreva a função em `lta/exportacao/` e registre-a em `FORMATOS` no `__init__.py`. Ela aparece automaticamente na aba Exportar.
- **Nova estrutura geométrica:** crie um construtor em `lta/geometria/construtores.py` que retorne um `Layout`. Todos os exportadores funcionam sem alteração.
