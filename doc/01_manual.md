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

## 3. Exportação

### 3.1 Layout

| Arquivo | Ferramenta | Como abrir |
|---|---|---|
| `.kicad_pcb` | KiCad 6/7/8 | Abrir no editor de PCB e pressionar **B** para preencher a zona de terra |
| `.kicad_pcb` | Altium | *File › Import Wizard › KiCad* (versões recentes; confirmar na versão instalada) |
| `.dxf` (R12) | Altium, ADS, HFSS, CST | Altium: *File › Import › DXF/DWG* no PCB, mapeando as camadas |
| `_topo.gtl`, `_terra.gbl`, `_contorno.gko` | Altium (CAMtastic), fabricação | Importar como Gerber RS-274X, em mm |
| `.gds` | ADS, HFSS, CST | ADS: *Layout › File › Import › GDSII*. Unidade de banco de dados 1 nm |
| `_hfss.py` | Ansys HFSS / AEDT | *Tools › Run Script*. Cria o projeto com variáveis `h_sub` e `h_ar`, portas de onda, radiação, setup e varredura |
| `_cst.bas` | CST Studio Suite | *Home › Macros* (editar/executar macro VBA). Parâmetros `h`, `t`, `er` e `tand` ficam editáveis |
| `_LEIA-ME.txt` | — | Substrato, camadas e parâmetros de projeto da estrutura |

Camadas (DXF/GDS): `TOPO_COBRE` (GDS 1), `TERRA` (GDS 2), `TERRA_FUROS` (GDS 3, furos a subtrair do terra) e `PLACA` (GDS 10).

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
| Exemplo Fig. 2.34: L, C, Cls | 3,1 nH, 0,17 pF, 1,2 pF | 3,30 nH, 0,157 pF, 1,31 pF (\*\*) |

(\*) O livro obtém εeff por solver EM; aqui ele vem de fórmula fechada (εeff 6,82 contra ≈ 7,0).
(\*\*) Os valores do livro já incluem otimização e arredondamento.

## 5. Limitações conhecidas

- Microstrip: Z0 e εeff com t = 0. Perdas e dispersão são opcionais (ver seção 2) e não incluem radiação, ondas de superfície nem a perda dos patches.
- Os scripts HFSS e CST ainda modelam o cobre como condutor perfeito (PEC); σ e rugosidade não são exportados.
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
