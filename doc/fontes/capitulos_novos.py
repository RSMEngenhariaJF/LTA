"""Seções e capítulos do manual acrescentados a partir da v1.2 do LTA (caps. 0 e 2 a 6 do livro).

Chamado por gerar_manual.py. Os números citados vêm dos testes automatizados (testes/test_validacao.py)
e dos exemplos de cada aba, reproduzidos nas capturas de tela.
"""
from pathlib import Path

IMG = Path(__file__).resolve().parent / "img"
LIVRO = "Martín (2015)"


def ref(secao, paginas=None):
    return f"{LIVRO}, seção {secao}" + (f", pp. {paginas}" if paginas else "")


# =============================================================================
# Acréscimos aos capítulos 2, 3 e 6 (linhas convencionais, Floquet, aplicações EBG)
# =============================================================================
def cap2_extras(m):
    m.h("2.3  Linhas planares convencionais: a ferramenta do capítulo 0", 2)
    m.p("A microstrip não é a única linha usada em placas de circuito impresso. Como as linhas artificiais "
        "dos capítulos seguintes são construídas sobre uma linha hospedeira, o LTA reúne no grupo "
        "**Cap. 0 – Linhas convencionais** o projeto das linhas planares mais comuns:")
    m.tabela(["Linha", "Geometria", "Modelo usado no LTA"],
             [["Microstrip", "trilha W sobre substrato h, terra embaixo",
               "Hammerstad e Jensen (1980), com dispersão e perdas opcionais"],
              ["Stripline", "fita W centrada entre dois planos de terra separados por b",
               "Cohn (1954), exata por transformação conforme (t = 0)"],
              ["CPW (coplanar)", "fita central S com fendas G e planos de terra no mesmo lado",
               "Simons (2001), integrais elípticas, substrato de espessura h"],
              ["CPW com terra (GCPW)", "CPW com plano de terra também na face inferior",
               "Simons (2001); cerca de vias opcional no layout"],
              ["CPS (coplanar strips)", "duas fitas paralelas separadas por S, sem terra",
               "Simons (2001); dual da CPW"]],
             [3.4, 6.2, 6.4], "2.1", "Linhas planares disponíveis na aba Cap. 0.")
    m.p("Nas linhas coplanares, a impedância depende da razão k = S/(S + 2G) por meio da integral elíptica "
        "completa K. Para a CPW, com k' = √(1 − k²):")
    m.eq(r"Z_{0}=\frac{30\pi}{\sqrt{\varepsilon_{\mathrm{eff}}}}\frac{K(k')}{K(k)},\qquad"
         r"\varepsilon_{\mathrm{eff}}=1+\frac{\varepsilon_{r}-1}{2}\frac{K(k_{1})K(k')}{K(k_{1}')K(k)}", "2.10")
    m.p("em que k₁ leva em conta a espessura finita h do substrato. Quando h cresce muito, εeff tende a "
        "(εr + 1)/2, valor que o LTA reproduz com erro de 10⁻⁵. Os testes também conferem a dualidade de "
        "Booker entre CPW e CPS, Z_CPW·Z_CPS = η₀²/(4εeff), com erro de 0,14%.")
    m.p("Além do trecho de linha com comprimento elétrico escolhido, a aba projeta, em microstrip, três "
        "elementos clássicos:")
    m.lista(["**Stub em aberto** em derivação, Y_in = j·tan(βl)/Z₀: com 90° em f₀ produz um zero de "
             "transmissão (rejeita-faixa).",
             "**Stub em curto** em derivação, Y_in = −j/(Z₀·tan βl): com 90° em f₀ é transparente em f₀ "
             "e bloqueia DC.",
             "**Transformador λ/4** entre Z₀ e uma carga Z₂, com Z_t = √(Z₀·Z₂). O cálculo usa portas de "
             "impedâncias diferentes, e o texto informa a banda em que |S₁₁| < −20 dB."])
    m.figura(IMG / "tela_17_linhas_gcpw_estrutura.png", "2.2",
             "Aba Linhas de transmissão convencionais: GCPW de 50 Ω em RO4003C (h = 0,813 mm) com cerca de "
             "vias. O corte A mostra a fita central, os planos coplanares, as vias e o terra inferior.")
    m.figura(IMG / "tela_18_linhas_transformador_resultados.png", "2.3",
             "Transformador λ/4 em microstrip entre 50 e 100 Ω a 2,45 GHz: |S₁₁| e |S₂₁| com portas de "
             "50 e 100 Ω e fase de transmissão.")
    m.p("CPW e CPS não têm plano de terra inferior. Nesse caso os arquivos 2D omitem a camada de terra, e "
        "nos scripts de HFSS e CST o terra fica inerte; as portas de onda precisam ser ajustadas, como "
        "indicado no arquivo LEIA-ME gerado. A stripline, por ter três camadas, só tem análise.")

    m.h("2.4  Perdas e dispersão", 2)
    m.p("As fórmulas da Seção 2.1 são quase-estáticas e sem perdas. Para estimar a perda de inserção real e "
        "o deslocamento das bandas em frequência alta, as abas EBG e Carga capacitiva têm a caixa **Modelo "
        "físico da linha**, com duas opções independentes (" + ref("1.5") + "):")
    m.lista(["**Perdas:** a constante de propagação passa a ser γ = α + jβ, com α = α_c + α_d. A perda "
             "dielétrica segue a Eq. (2.11), e a do condutor usa a resistência superficial R_s = 1/(σδ), "
             "com correções de espessura finita e de rugosidade (Hammerstad).",
             "**Dispersão:** εeff(f) de Kirschning e Jansen (1982), válido para εr ≤ 20 e f·h ≤ 25 GHz·mm, "
             "e Z₀(f) de Bianco et al. (1978)."])
    m.eq(r"\alpha_{d}=\frac{k_{0}\,\varepsilon_{r}(\varepsilon_{\mathrm{eff}}-1)\tan\delta}"
         r"{2\sqrt{\varepsilon_{\mathrm{eff}}}(\varepsilon_{r}-1)},\qquad"
         r"\delta=\frac{1}{\sqrt{\pi f\mu_{0}\sigma}}", "2.11")
    m.eq(r"K_{\mathrm{rug}}=1+\frac{2}{\pi}\arctan\left[1,4\left(\frac{\Delta}{\delta}\right)^{2}\right]", "2.12")
    m.p("A condutividade σ e a rugosidade Δ são definidas na aba Substrato. Como referência, a profundidade "
        "pelicular do cobre calculada pelo LTA é 2,09 µm a 1 GHz e 0,66 µm a 10 GHz, contra 2,06 e 0,66 µm "
        "na nota 24 do capítulo 1 do livro.")
    m.figura(IMG / "tela_37_ebg_perdas.png", "2.4",
             "EBG senoidal da Fig. 2.11 com perdas e dispersão. O centro da banda passa de 5,82 para "
             "5,64 GHz (livro, com solver EM: 5,66 GHz). O último gráfico mostra α(f) e εeff(f) da linha de "
             "referência.")
    m.p("Com a dispersão ligada, o Z₀ das pontas da linha também cresce com a frequência (cerca de 54 Ω a "
        "10 GHz em RO3010). Por isso as curvas de |S₁₁| da cascata ABCD, referidas a 50 Ω fixos, e do CMT, "
        "referidas ao Z₀ local, se afastam em frequência alta. As soluções analíticas do Capítulo 6 "
        "continuam sem perdas e sem dispersão.")


def cap3_extras(m):
    m.h("3.3  Harmônicos espaciais na prática: a ferramenta Floquet", 2)
    m.p("A aba **Floquet / harmônicos espaciais** desenha o diagrama de dispersão na zona estendida, como na "
        "Fig. 2.2 do livro. O βl calculado pela matriz ABCD fica na zona reduzida [0, π]; a ferramenta o "
        "**desdobra**, escolhendo entre 2πk ± βl o valor contínuo e não decrescente. Cada harmônico espacial "
        "tem constante de fase deslocada de 2π/l (" + ref("2.2") + "):")
    m.eq(r"\beta_{n}l=\beta l+2\pi n,\qquad n=0,\pm1,\pm2,\ldots", "3.4")
    m.p("O harmônico n radia quando é uma onda rápida, |β_n| < k₀. O ângulo do feixe, medido a partir do "
        "eixo da linha (90° = *broadside*), é dado por:")
    m.eq(r"\cos\theta_{n}=\frac{\beta_{n}}{k_{0}}", "3.5")
    m.p("A ferramenta sombreia a região de onda rápida, destaca os trechos que radiam e lista as faixas de "
        "radiação de cada harmônico: para trás (β_n < 0) ou para a frente (β_n > 0). Ela também estima onde "
        "começa a radiação do harmônico n = −1 e onde o feixe passa por *broadside*. Os exemplos são a EBG da "
        "Fig. 2.11, a linha uniforme (retas deslocadas de 2π/l, Fig. 2.1) e a linha com carga capacitiva.")
    m.figura(IMG / "tela_19_floquet.png", "3.2",
             "Aba Floquet com a EBG da Fig. 2.11: diagrama na zona estendida com a região de onda rápida, "
             "atenuação por célula, ângulo do feixe (Eq. 3.5) e velocidades de fase e de grupo.")
    m.nota("a radiação não entra no cálculo de α. A ferramenta indica **onde** a estrutura pode radiar, "
           "mas não **quanto** se perde; para isso é preciso simulação eletromagnética. O desdobramento supõe "
           "estrutura de mão direita; as linhas CRLH do Capítulo 7 usam o sinal de β calculado na própria "
           "aba CRLH.", "Limite do modelo")


def cap6_extras(m):
    m.h("6.8  Outras aplicações EBG: linha de atraso chirpada e ressonador de Bragg", 2)
    m.p("A aba **Aplicações EBG** projeta duas estruturas da seção 2.5.1 do livro a partir de "
        "especificações, usando o mesmo CMT numérico da aba EBG.")
    m.h("Linha de atraso chirpada (CDL)", 3)
    m.p("Numa linha chirpada, o período da modulação varia ao longo de z. Cada frequência é refletida numa "
        "posição diferente, e o atraso de grupo da reflexão varia linearmente com a frequência, com "
        "inclinação ψ (em ns/GHz). Esse comportamento é útil em compressão de pulsos e em análise de "
        "espectro em tempo real (" + ref("2.5.1.5") + "). O perfil de impedância é:")
    m.eq(r"Z_{0}(z)=Z_{\mathrm{ref}}\exp\left[A\,W(z)\sin\left(\frac{2\pi z}{l_{0}}+Cz^{2}-\frac{CL^{2}}{4}\right)\right]", "6.12")
    m.p("Dadas f₀, a banda Δf e ψ, o LTA calcula o período central l₀, o coeficiente de chirp C e o "
        "comprimento L (Eqs. 2.109–2.115 do livro). L é arredondado para um número inteiro de períodos, para "
        "que as portas fiquem em Z_ref, e C é recalculado para manter Δf. A janela W(z) pode ser a gaussiana "
        "assimétrica do livro, Hamming ou nenhuma.")
    m.p("No exemplo do livro (Fig. 2.26), o LTA obtém l₀ = 6,39 mm e C = −2052 m⁻² (livro: 6,4 mm e "
        "−2080 m⁻²) e uma inclinação simulada sem perdas de −0,48 ns/GHz (projeto: −0,46). Com perdas, a "
        "perda média de reflexão é de cerca de 3,4 dB (livro: cerca de 3 dB).")
    m.figura(IMG / "tela_20_cdl.png", "6.8",
             "Linha de atraso chirpada da Fig. 2.26 (N = 50), com perdas e dispersão ligadas: perfil Z₀(z) "
             "e largura, frequência refletida em cada posição, |S₁₁| e |S₂₁| e atraso de grupo de S₁₁ "
             "comparado com a reta de projeto.")
    m.p("As Eqs. 2.111–2.115 supõem εeff constante. Com a dispersão ligada, a banda real se comprime e ψ "
        "muda (cerca de +19% no exemplo), o que deve ser compensado no projeto.")
    m.h("Ressonador de Bragg", 3)
    m.p("Dois refletores EBG senoidais sintonizados em f₀, separados por uma cavidade uniforme de "
        "comprimento L_m, formam um ressonador de alto Q (" + ref("2.5.1.2") + "). As fases de "
        "reflexão dos refletores, φ₁ e φ₂, vêm do CMT em f₀, e L_m sai da condição de ressonância:")
    m.eq(r"\varphi_{1}+\varphi_{2}-2\beta L_{m}=2\pi n", "6.13")
    m.p("A ferramenta escolhe o n que dá o menor L_m positivo e calcula o pico de |S₂₁|, a perda de inserção, "
        "o Q carregado Q_L = f₀/BW₋₃dB e o Q descarregado Q₀ = Q_L/(1 − 10^(−IL/20)). No exemplo da Fig. 2.21 "
        "(4,2 GHz em RO3010), o pico cai exatamente em f₀ e Q_L = 132 com perdas. O livro mede Q = 129,2, mas "
        "com refletores de furos no terra, que não têm modelo de circuito no LTA.")
    m.figura(IMG / "tela_21_ressonador.png", "6.9",
             "Ressonador de Bragg da Fig. 2.21: perfil de impedância dos dois refletores e da cavidade, "
             "resposta pelo CMT, detalhe da ressonância (Q_L = 132,3) e fases de reflexão dos refletores.")
    m.p("O filtro *wiggly-line* (seção 2.5.1.3) não está implementado, porque exige um modelo de linhas "
        "acopladas.")


# =============================================================================
# Capítulo 7: metamateriais (cap. 3 do livro)
# =============================================================================
def cap7_metamateriais(m):
    m.h("7  Linhas metamateriais: CRLH, CSRR e extração de parâmetros")
    m.p("O capítulo 3 do livro troca o ponto de vista: em vez de perturbar uma linha comum, a célula é "
        "projetada como um circuito LC que pode ter comportamento de **mão esquerda** (*left-handed*, LH), "
        "com velocidades de fase e de grupo em sentidos opostos. Isso acontece quando o ramo série é "
        "capacitivo e o ramo em derivação é indutivo. Combinada à linha hospedeira, que é de mão direita "
        "(RH), essa célula forma a linha **CRLH** (*composite right/left-handed*) (" + ref("3.4", "149–162") + ").")
    m.p("No capítulo 3 não existe fórmula fechada confiável que leve do circuito à geometria. O livro fecha "
        "os projetos com simulação eletromagnética e extração de parâmetros (Apêndice G) ou com *space "
        "mapping* (Apêndice H). As quatro ferramentas do grupo **Cap. 3 – Metamateriais** seguem esse fluxo:")
    m.lista(["**circuito:** síntese dos elementos a partir das especificações (aba CRLH);",
             "**geometria inicial:** dimensões estimadas e layout (abas CRLH CL-loaded e Linha CSRR + gap);",
             "**simulação EM ou medição** da célula exportada (HFSS, CST);",
             "**extração** dos elementos a partir do .s2p simulado ou medido (aba Extração);",
             "**correção** da geometria com o fator de calibração obtido, e nova iteração."],
            numerada=True)

    m.h("7.1  A linha CRLH", 2)
    m.p("A célula CRLH tem, no ramo série, L_R em série com C_L e, em derivação, C_R em paralelo com L_L. "
        "Sua dispersão vem da mesma relação de Bloch do Capítulo 4, com as impedâncias Z_s do ramo série e "
        "Y_p do ramo em derivação:")
    m.eq(r"\cos\beta l=1+\frac{Z_{s}Y_{p}}{2},\qquad Z_{s}=j\left(\omega L_{R}-\frac{1}{\omega C_{L}}\right),"
         r"\quad Y_{p}=j\left(\omega C_{R}-\frac{1}{\omega L_{L}}\right)", "7.1")
    m.p("Em baixa frequência dominam C_L e L_L e a linha é LH (β < 0); em alta frequência dominam L_R e C_R e "
        "ela é RH. Entre as duas regiões há, em geral, uma banda proibida limitada pelas ressonâncias série "
        "ω_s = 1/√(L_R·C_L) e em derivação ω_p = 1/√(L_L·C_R). Quando ω_s = ω_p, a linha é **balanceada**: "
        "a banda proibida desaparece e a transição LH → RH ocorre em ω₀ com β = 0 e impedância de Bloch "
        "finita:")
    m.eq(r"L_{R}C_{L}=L_{L}C_{R},\qquad Z_{B}=\sqrt{\frac{L_{R}}{C_{R}}}=\sqrt{\frac{L_{L}}{C_{L}}},\qquad"
         r"\beta l\approx\frac{\omega}{\omega_{R}}-\frac{\omega_{L}}{\omega}", "7.2")
    m.p("com ω_R = 1/√(L_R·C_R) e ω_L = 1/√(L_L·C_L). A aba **CRLH / D-CRLH** sintetiza a célula "
        "balanceada de duas formas: pela banda passante (bordas f_CL e f_CR e Z_B, Eq. 3.66 do livro) ou "
        "pelas fases em duas frequências, que é o ponto de partida dos dispositivos dual-band do Capítulo 8. "
        "Também analisa células quaisquer e a D-CRLH (dual), que tem banda proibida no lugar da banda "
        "passante.")
    m.figura(IMG / "tela_22_crlh.png", "7.1",
             "Aba CRLH com a célula não balanceada da Fig. 3.21: dispersão com sinal (LH em vermelho, RH em "
             "azul) e região de onda rápida, impedância de Bloch, parâmetros S de 7 células, εeff e μeff "
             "efetivos, fase de transmissão e ângulo do feixe leaky-wave do modo fundamental.")
    m.p("Na faixa em que |β| < k₀, o modo fundamental da linha CRLH é uma onda rápida e pode radiar. O feixe "
        "varre de *backfire* (banda LH) a *endfire* (banda RH), passando por *broadside* na transição: é o "
        "princípio das antenas *leaky-wave* de varredura contínua, tratadas no Capítulo 8.")
    m.tabela(["Ao aumentar…", "Efeito na linha CRLH", "Impacto prático"],
             [["C_L (série)", "desce ω_s e a borda inferior da banda LH", "banda LH mais larga para baixo; "
               "capacitores interdigitais maiores"],
              ["L_L (derivação)", "desce ω_p", "stubs mais longos; desbalanceia se não acompanhado de C_L"],
              ["L_R e C_R (hospedeira)", "descem ω_R e a borda superior da banda RH",
               "célula mais longa; menos banda total"],
              ["Desbalanceamento ω_s ≠ ω_p", "abre uma banda proibida entre LH e RH",
               "Z_B vai a zero ou a infinito na borda; feixe não passa por broadside"]],
             [3.6, 5.6, 6.8], "7.1", "Influência dos elementos da célula CRLH.")

    m.h("7.2  Implementação em microstrip: CRLH CL-loaded", 2)
    m.p("A forma mais direta de construir uma linha CRLH em microstrip é carregar a linha hospedeira com "
        "capacitores série e indutores em derivação (" + ref("3.5.1") + ", Fig. 3.27). A aba **CRLH "
        "microstrip (CL-loaded)** converte os quatro elementos em geometria:")
    m.lista(["a linha hospedeira de impedância √(L_R/C_R) e período p fornece L_R;",
             "o **capacitor interdigital** fornece C_L, calculado pela fórmula de Bahl; nas pontas da linha "
             "usam-se capacitores de 2C_L, para que a célula fique simétrica;",
             "o **stub em curto** de impedância escolhida, com a indutância da **via** (Goldfarb e Pucel), "
             "fornece L_L."])
    m.p("Os avisos da aba indicam quando o interdigital não cabe no período ou quando o stub fica longo "
        "demais. Isso não é uma limitação do programa, mas a limitação real de banda larga em uma única "
        "camada, discutida na seção 3.5.1 do livro. Com L_R = 2,45 nH, C_R = 0,5 pF, L_L = 3,38 nH e "
        "C_L = 0,68 pF, a transição fica em 3,885 GHz (livro, Fig. 3.28: 3,9 GHz).")
    m.figura(IMG / "tela_23_cl_loaded_estrutura.png", "7.2",
             "Aba CRLH CL-loaded: linha hospedeira com capacitores interdigitais em série e stubs com via em "
             "derivação; o corte mostra a via metalizada.")

    m.h("7.3  Linha com CSRR e gaps", 2)
    m.p("A segunda família de linhas metamateriais usa **ressonadores**. Um CSRR (*complementary split-ring "
        "resonator*) gravado no plano de terra, sob a trilha, se comporta como um tanque L_c‖C_c acoplado à "
        "linha pela capacitância C. Um gap série na trilha fornece a capacitância C_g. O resultado é uma "
        "célula CRLH ressonante, cujo circuito equivalente é a Fig. 3.38b do livro "
        "(" + ref("3.5.2.2") + "). Para a célula balanceada com transição em f₀, impedância Z_B e zero de "
        "transmissão em f_z, o LTA usa:")
    m.eq(r"C_{c}=\frac{L}{Z_{B}^{2}},\quad L_{c}=\frac{1}{\omega_{0}^{2}C_{c}},\quad"
         r"C_{g}=\frac{1}{\omega_{0}^{2}L},\quad C=\frac{1}{\omega_{z}^{2}L_{c}}-C_{c}", "7.3")
    m.p("em que L é a indutância da linha no período escolhido. O **modelo melhorado** (Eqs. 3.87–3.88) "
        "separa a capacitância série física do gap, C_s = C_g²/(C + 2C_g), da parcela paralela "
        "C_par = C_g·C/(C + 2C_g). Nos valores da Tabela 3.3, C_s = 0,014 pF e C_par = 0,28 pF: o C_g do "
        "modelo é quase todo capacitância paralela, como o livro argumenta.")
    m.p("O raio do CSRR é estimado por um modelo quase-estático (seção 3.3.2.1 do livro). Na Fig. 3.39, essa "
        "estimativa dá a ressonância 16% abaixo da extraída. Por isso a aba aceita um **fator de "
        "calibração** (f_o real / f_o estimada), obtido de uma simulação ou medição pela aba Extração; com ele, "
        "o raio de 3,661 mm do livro é reproduzido.")
    m.figura(IMG / "tela_24_csrr_estrutura.png", "7.3",
             "Aba Linha CSRR + gap com a célula balanceada da Fig. 3.39: trilha com o gap interdigital no "
             "topo e o CSRR gravado no terra (tracejado na vista superior).")
    m.figura(IMG / "tela_25_csrr_resultados.png", "7.4",
             "Resultados da linha CSRR + gap: dispersão, impedância de Bloch e parâmetros S do circuito "
             "balanceado.")

    m.h("7.4  Extração de parâmetros a partir de simulação ou medição", 2)
    m.p("A extração (Apêndice G do livro) é o elo entre o circuito e a geometria. A aba **Extração de "
        "parâmetros** lê um arquivo .s2p, simulado em EM ou medido, e determina os elementos do circuito "
        "equivalente usando condições que dependem só da topologia:")
    m.lista(["S₁₁ cruza o círculo **r = 1** da carta de Smith onde o ramo em derivação abre;",
             "S₁₁ cruza o círculo **g = 1** onde a reatância do ramo série se anula;",
             "o **zero de transmissão** dá a ressonância do ramo em derivação;",
             "o **zero de reflexão** e a frequência em que βl = ±π/2 completam o sistema."])
    m.p("Como os cruzamentos podem ser tangências ou estar fora da faixa, o LTA reúne todos os candidatos "
        "(cruzamentos e zeros dos equivalentes T e π da célula medida), testa as combinações e fica com a de "
        "menor erro rms. Erro abaixo de 0,02 indica bom ajuste; acima de 0,1, modelo inadequado ou faixa "
        "grande demais. Estão disponíveis sete modelos: CSRR, CSRR + gap, SRR, SRR + fita, OSRR, OCSRR e "
        "CRLH CL-loaded. Na ida e volta com os valores das Tabelas 3.1 e 3.3, os elementos são recuperados "
        "com desvio menor que 0,5% (menor que 2% com perdas no ressonador).")
    m.figura(IMG / "tela_26_extracao.png", "7.5",
             "Aba Extração com o exemplo CSRR + gap da Tabela 3.3: S₁₁ na carta de Smith com os círculos "
             "r = 1 e g = 1, módulo e fase de S e dispersão da célula. Dados e modelo extraído se "
             "sobrepõem; os elementos e as frequências usadas aparecem no painel.")
    m.nota("quando a linha CSRR + gap é balanceada (f_s ≈ f_o), a reatância série se anula em f_o e não "
           "determina L. O LTA usa então Z_B(βl = π/2) = |X do braço série|, que é exato para a célula "
           "simétrica sem perdas. Na transformação SRR + fita das Eqs. 3.71–3.74, o fator (1 + M²/2L_pL_s) "
           "só reproduz a Tabela 3.1 se estiver elevado ao quadrado; o expoente não aparece no texto extraído "
           "do PDF.", "Detalhes de implementação")


# =============================================================================
# Capítulo 8: aplicações (cap. 4 do livro)
# =============================================================================
def cap8_aplicacoes(m):
    m.h("8  Aplicações: dispositivos multibanda, híbridas, filtros e antenas")
    m.p("O capítulo 4 do livro usa as linhas artificiais para resolver um problema que a linha comum não "
        "resolve: obter **a fase certa em mais de uma frequência**. As ferramentas do grupo **Cap. 4 – "
        "Aplicações** trabalham no nível de circuito; o caminho até o layout é o mesmo do Capítulo 7 "
        "(" + ref("4.2", "215–290") + ").")

    m.h("8.1  Linhas dual-band e multibanda", 2)
    m.p("Uma linha convencional de 90° em f₁ tem 270° em 3f₁, e só nessas frequências harmônicas um "
        "dispositivo λ/4 funciona de novo. Uma célula artificial tem quatro elementos livres e pode impor "
        "fase e impedância em duas frequências quaisquer. Para uma célula T com fase φ e impedância Z, os "
        "braços série e o ramo em derivação valem:")
    m.eq(r"X_{a}=Z\tan\frac{\varphi}{2},\qquad X_{p}=-\frac{Z}{\sin\varphi}", "8.1")
    m.p("Impondo, por exemplo, φ = −90° em f₁ e +90° em f₂, saem quatro equações para os elementos. A aba "
        "**Linhas dual-band e multibanda** resolve quatro topologias:")
    m.lista(["**linha convencional + PLH** (linha de mão esquerda pura), Eqs. 4.21–4.24. As fórmulas supõem "
             "a PLH homogênea: com uma célula a fase erra 23° em f₁; com quatro, 2°;",
             "**CRLH canônica** (Eqs. 4.29), balanceada em √(ω₁ω₂) quando Z₁ = Z₂;",
             "**CSRR + gap**, com zero de transmissão e trechos de linha nas pontas (Eqs. 4.33–4.34). Os dois "
             "exemplos do livro são reproduzidos com três algarismos (Tabela 8.1);",
             "**E-CRLH quad-band** (Eqs. 4.37): oito condições, ±90° e 35,35 Ω em quatro frequências."])
    m.tabela(["Caso (livro)", "Grandeza", "Livro", "LTA"],
             [["CSRR + gap, p. 258", "L, C, C_g, L_c, C_c", "12,5 / 24,9 / 1,25 / 3,38 / 5,10",
               "12,50 / 24,91 / 1,251 / 3,376 / 5,096"],
              ["Idem com trechos de linha, φ_L(f₁) = 15°", "L, C, C_g, L_c, C_c",
               "9,45 / 17,9 / 1,01 / 4,85 / 2,95", "9,45 / 17,92 / 1,005 / 4,852 / 2,956"],
              ["Inversor CRLH 2,4/3,75 GHz, p. 255", "C_L, L_L", "0,33 pF; 0,83 nH", "0,338 pF; 0,844 nH"],
              ["Quad-band E-CRLH, p. 267", "8 elementos", "valores publicados",
               "a até 3,4% (fases e Z exatas)"]],
             [4.4, 3.0, 4.3, 4.3], "8.1", "Projetos dual-band e quad-band: comparação com o livro "
             "(indutâncias em nH, capacitâncias em pF).", fonte=8.5)
    m.p("No caso quad-band, a solução do LTA satisfaz exatamente as oito condições, enquanto os valores "
        "publicados erram a impedância em até 4 Ω; a topologia da Fig. 3.25c foi identificada a partir "
        "desses valores.")
    m.figura(IMG / "tela_27_multibanda_quad.png", "8.1",
             "Aba Linhas dual-band e multibanda com a célula E-CRLH quad-band (0,9; 1,17; 1,57 e 1,8 GHz): "
             "fase de transmissão com os alvos de ±90° (pontos), impedância de Bloch e S com portas de "
             "35,35 Ω.")

    m.h("8.2  Híbridas e divisores", 2)
    m.p("A aba **Híbridas e divisores** monta branch-line, rat-race, Wilkinson e divisor com inversor como "
        "redes de N portas: cada linha entra na matriz Y nodal pela sua matriz ABCD, e os nós internos são "
        "eliminados por redução de Kron. As linhas podem ser convencionais, CRLH dual-band ou, no rat-race, "
        "uma linha PLH de −90° no lugar do trecho de 270° (Figs. 4.1 e 4.5 do livro). Como a PLH tem fase "
        "que varia em sentido oposto ao da linha comum, a diferença de 180° entre os caminhos se mantém numa "
        "faixa muito maior.")
    m.p("A ferramenta mede a faixa contínua em torno de cada frequência de projeto em que valem os critérios "
        "da Tabela 4.1: |S₁₁| < −10 dB, saídas em −3 ± 1 dB, balanço de fase ± 5° e isolamento < −15 dB. No "
        "rat-race de banda larga, o balanço de fase do circuito ideal passa de 15,6% para 102% de banda.")
    m.figura(IMG / "tela_28_ratrace.png", "8.2",
             "Rat-race de banda larga com linha LH de −90° (esquerda) contra o convencional em 1,6 GHz "
             "(direita): parâmetros S, balanço de fase e potência entregue às saídas.")

    m.h("8.3  Filtros passa-faixa e dual-band", 2)
    m.p("A aba **Filtros** parte dos protótipos passa-baixas de Butterworth ou Chebyshev (coeficientes g) e "
        "projeta filtros passa-faixa com inversores de impedância e ressonadores em derivação "
        "(" + ref("4.2.3") + "). Cada tanque tem:")
    m.eq(r"C_{eq}=\frac{g_{i}}{\mathrm{FBW}\,\omega_{0}Z_{0}},\qquad L_{eq}=\frac{1}{\omega_{0}^{2}C_{eq}}", "8.2")
    m.p("A **transformação dual-band** (Eqs. 4.49–4.50) troca cada tanque por (L'‖C')‖(L'' + C''), que é a "
        "derivação da célula E-CRLH, e produz duas bandas passantes. A terceira opção é o **filtro periódico "
        "com células híbridas de CSRR**, que usa a mesma célula repetida: C_g pela Eq. 4.46 e os demais "
        "elementos pelas Eqs. 3.86 e 4.43–4.45. Com três células e elementos iguais, g = 1,521, como no livro.")
    m.figura(IMG / "tela_29_filtro_hibrido.png", "8.3",
             "Filtro periódico de três células híbridas de CSRR a 1 GHz: resposta do circuito com FBW de "
             "9,5% (o protótipo medido do livro tem 8%).")
    m.nota("os valores da Tabela 4.2 do livro dão zero de transmissão em 1,80 GHz e uma reatância em "
           "derivação de j33,7 Ω em f₀. Provavelmente já incluem o ajuste de layout descrito no texto, por "
           "isso a diferença em relação à síntese do LTA é esperada.", "Diferença com o livro")

    m.h("8.4  Antenas leaky-wave", 2)
    m.p("Na aba CRLH, informando a constante de atenuação normalizada α/k₀ (obtida em simulação EM), o LTA "
        "calcula o número de células necessário para irradiar uma fração η da potência e estima a largura do "
        "feixe pela aproximação de abertura de Jackson e Oliner:")
    m.eq(r"\eta=1-e^{-2\alpha L}\;\Rightarrow\;L=\frac{-\ln(1-\eta)}{2\alpha}", "8.3")
    m.p("O ângulo do feixe vem de cos θ = β/k₀ (Eq. 3.5), com o β do modo fundamental. Numa linha CRLH "
        "balanceada, ele varre continuamente de *backfire* a *endfire* ao variar a frequência; no Capítulo 9, "
        "a mesma varredura é obtida em frequência fixa, sintonizando a linha.")


# =============================================================================
# Capítulo 9: sintonizáveis e não lineares (cap. 5 do livro)
# =============================================================================
def cap9_sintonia(m):
    m.h("9  Linhas sintonizáveis e não lineares")
    m.p("Quando um dos elementos da célula pode ser controlado por tensão, a linha artificial passa a ser "
        "**sintonizável**: suas bandas, sua fase ou seu feixe de radiação mudam sem alterar a geometria. "
        "Quando a capacitância depende da própria tensão do sinal, a linha é **não linear** e pode formar "
        "sólitons, comprimir pulsos e gerar harmônicos (" + LIVRO + ", capítulo 5).")

    m.h("9.1  Elementos de sintonia", 2)
    m.p("A aba **Linhas sintonizáveis** traz os modelos dos três dispositivos tratados no livro:")
    m.eq(r"C_{\mathrm{var}}(V)=\frac{C_{j0}}{(1+V/\phi)^{\gamma}},\qquad C_{\mathrm{HBV}}(V)=C_{M}\,"
         r"\mathrm{sech}(kV),\qquad V_{p}=\sqrt{\frac{8kg_{0}^{3}}{27\varepsilon_{0}A}}", "9.1")
    m.lista(["**varactor:** junção com γ entre 1/3 e 1/2;",
             "**HBV** (*heterostructure barrier varactor*): capacitância par em V, útil para gerar harmônicos "
             "ímpares;",
             "**RF-MEMS** de placas paralelas: C = ε₀A/g; a placa colapsa (*pull-in*) na tensão V_p, quando "
             "o vão chega a 2g₀/3 (Eq. 5.5 do livro). O LTA confere g = 0,6667·g₀ no colapso."])

    m.h("9.2  Rejeita-faixa sintonizável com SRR e varactor", 2)
    m.p("Um SRR acoplado magneticamente à linha, com um varactor no lugar da capacitância do anel (VLSRR), "
        "produz um tanque em série refletido na linha (Eqs. 5.6). Com C_eq = C_L em série com C_var:")
    m.eq(r"L_{s}'=\frac{4M^{2}}{L_{r}},\qquad C_{s}'=\frac{L_{r}}{4M^{2}\omega_{0}^{2}},\qquad"
         r"R=\frac{4M^{2}\omega_{0}^{2}}{R_{s}},\qquad Q=\frac{\omega_{0}L_{r}}{R_{s}}", "9.2")
    m.p("Com os valores da Fig. 5.7, a rejeição se desloca de 2,64 a 3,80 GHz e fica mais funda no topo da "
        "faixa (de −32 para −48 dB), como o livro descreve.")
    m.figura(IMG / "tela_30_vlsrr.png", "9.1",
             "Linha com SRR e varactor (Fig. 5.7): |S₂₁| de duas células para C_var de 1 a 9 pF e frequência "
             "de rejeição em função de C_var.")

    m.h("9.3  CRLH sintonizável: defasador e antena de varredura", 2)
    m.p("Se C_L e L_L forem multiplicados pelo mesmo fator x (por exemplo, com varactores e indutores "
        "sintonizáveis), a linha CRLH continua balanceada, com Z_B constante, e a transição se desloca para "
        "f₀/√x (Eqs. 5.9–5.11). O LTA confirma Z_B = 50,00 Ω de x = 0,6 a 1,7. Duas aplicações decorrem:")
    m.lista(["**defasador:** a defasagem obtida numa frequência de operação fixa, com a figura de mérito da "
             "Eq. 5.12 (excursão de fase em graus dividida pela maior perda, em dB), considerando Q finito nos elementos sintonizáveis;",
             "**antena leaky-wave de varredura em frequência fixa:** o ângulo do feixe muda com a tensão, e não "
             "com a frequência (seção 5.3.2.2)."])
    m.figura(IMG / "tela_31_defasador.png", "9.2",
             "CRLH sintonizável como defasador: dispersão e |S₂₁| de seis células (Q = 40) para cinco "
             "fatores de sintonia, Re Z_B constante na transição e defasagem em 2,45 GHz em função de x.")

    m.h("9.4  Linhas de transmissão não lineares (NLTL)", 2)
    m.p("Numa linha LC periódica carregada com capacitâncias que dependem da tensão, os pulsos de amplitude "
        "maior viajam mais depressa (a capacitância cai com a tensão). A dispersão da estrutura periódica "
        "compensa esse efeito e surgem os **sólitons**: pulsos que se propagam sem mudar de forma, com atraso "
        "por seção que depende da amplitude (" + ref("5.4") + "). O LTA calcula:")
    m.lista(["a faixa de atrasos possíveis, de √(L·C(∞)) a √(L·C(0)): 2,174 a 4,216 ps na linha da Fig. 5.40;",
             "o atraso em função da amplitude (Eq. 5.24) e a forma do sóliton, por quadratura da Eq. 5.20;",
             "a **simulação no tempo** das Eqs. 5.13–5.15 (Runge-Kutta de 4ª ordem) com fonte e carga "
             "resistivas, para várias linhas de comprimentos diferentes de uma vez."])
    m.p("O atraso dos sólitons simulados confere com a Eq. 5.24 com erro de até 1,5%; o erro cresce para "
        "sólitons curtos, como esperado da aproximação contínua. A simulação também mostra a **decomposição "
        "de um pulso** largo em vários sólitons e a **multiplicação de frequência**: com os parâmetros da "
        "Fig. 5.44 e C_M = 0,12 pF, o 3º harmônico é máximo com 33 seções (livro: 27). A posição exata "
        "depende das terminações, que o livro não informa; o padrão no LTA é √(L/C₀) = 50 Ω.")
    m.figura(IMG / "tela_32_nltl_pulso.png", "9.3",
             "Aba NLTL: pulso de 10 V e 40 ps se decompondo em sólitons ao longo de 95 seções com HBVs "
             "(cf. Fig. 5.42b do livro).")


# =============================================================================
# Capítulo 10: outras linhas avançadas (cap. 6 do livro)
# =============================================================================
def cap10_avancadas(m):
    m.h("10  Outras linhas avançadas")
    m.p("O capítulo 6 do livro reúne linhas artificiais com aplicações específicas. O LTA implementa quatro "
        "delas no grupo **Cap. 6 – Outras linhas avançadas**.")

    m.h("10.1  Linhas diferenciais com supressão do modo comum", 2)
    m.p("Em interconexões diferenciais de alta velocidade, o ruído de modo comum causa interferência "
        "eletromagnética. CSRRs gravados no terra sob o par diferencial resolvem o problema de forma "
        "seletiva (" + ref("6.3") + "): no **modo diferencial** o plano de simetria é uma parede elétrica, o "
        "CSRR não é excitado e o sinal passa; no **modo comum** a parede é magnética, o CSRR ressoa e "
        "rejeita uma banda. O zero de transmissão do modo comum é:")
    m.eq(r"f_{z}=\frac{1}{2\pi\sqrt{2L_{c}\left(C_{e}+C_{c}/2\right)}}", "10.1")
    m.p("A aba calcula a **banda máxima de rejeição**, que é a faixa em que nenhum dos dois modos de Bloch da "
        "célula de quatro portas se propaga (Eqs. 3.96–3.97 do livro), e simula o filtro de N células. A banda "
        "a −20 dB converge para a máxima com o número de células; seis células já ficam perto dela, como o "
        "livro afirma.")
    m.tabela(["Linha (Tabela 6.1)", "f_z livro [GHz]", "f_z LTA [GHz]", "Banda máx. livro", "Banda máx. LTA"],
             [["CSRR 1,5 GHz", "1,5", "1,495", "30,7%", "32,7%"],
              ["CSRR 3 GHz", "3", "2,906", "26,4%", "26,9%"],
              ["CSRR 6 GHz", "6", "5,931", "20,7%", "21,6%"],
              ["DS-CSRR (três casos)", "–", "–", "32,6 / 29,4 / 23,6%", "32,4 / 28,9 / 21,9%"]],
             [4.0, 2.8, 2.8, 3.2, 3.2], "10.1", "Linha diferencial com CSRR: comparação com a Tabela 6.1 do "
             "livro (diferença máxima de 2,0 pontos percentuais na banda).")
    m.figura(IMG / "tela_33_diferencial.png", "10.1",
             "Aba Linha diferencial com o DS-CSRR de 3 GHz: S dos modos comum e diferencial para 6 células, "
             "atenuação de Bloch do modo comum e banda de rejeição em função do número de células.")

    m.h("10.2  Linhas de atraso com ondas magnetoindutivas (MIW/EIW)", 2)
    m.p("Uma cadeia de SRRs acoplados magneticamente (MIW) ou de CSRRs acoplados eletricamente (EIW) "
        "suporta ondas numa banda estreita em torno da ressonância dos anéis, com velocidade de grupo muito "
        "baixa (" + ref("6.2") + "). A dispersão é:")
    m.eq(r"\omega=\frac{\omega_{0}}{\sqrt{1+2\kappa\cos\beta l}},\qquad\kappa=\frac{M}{L}\;(\mathrm{MIW}),"
         r"\quad\kappa=-\frac{C_{M}}{C_{c}}\;(\mathrm{EIW})", "10.2")
    m.p("No centro da banda, v_g = l·κ·ω₀ e o atraso de ida e volta por ressonador é 2/(ω₀|κ|). A aplicação "
        "do livro é a **etiqueta RFID sem chip por reflectometria**: refletores ao longo da cadeia geram ecos "
        "cujos atrasos codificam bits. O LTA calcula o número de bits Δτ/τ e mostra o envelope dos ecos de "
        "um pulso gaussiano. Com |M|/L = 0,0256, a banda é de 5,13%, igual à do livro, e o atraso por "
        "ressonador confere com a cadeia simulada com erro menor que 0,1%.")
    m.figura(IMG / "tela_34_miw.png", "10.2",
             "Aba Linhas de atraso MIW/EIW com a etiqueta de 2,45 GHz do livro: dispersão (Eq. 6.1), fator "
             "de onda lenta c/|v_g| e resposta da etiqueta \"11\" a um pulso de 15 ns.")
    m.p("Pulsos de banda larga em relação à banda MIW chegam atrasados e geram ecos espúrios. Esse é o "
        "compromisso entre número de bits, banda, comprimento e perdas descrito no livro (p. 429).")

    m.h("10.3  Guia de onda integrado ao substrato (SIW)", 2)
    m.p("O SIW é um guia retangular construído numa placa de circuito impresso: duas fileiras de vias "
        "metalizadas fazem as paredes laterais (" + ref("6.5") + "). Ele herda o baixo nível de perdas e a "
        "blindagem do guia de onda, com a fabricação barata da placa. O modo TE₁₀ tem corte em f_c, e a "
        "largura efetiva e a distância entre fileiras valem:")
    m.eq(r"a=\frac{c}{2f_{c}\sqrt{\varepsilon_{r}}},\qquad W=\frac{1}{2}\left[a+\sqrt{(a+0,54d)^{2}-0,4d^{2}}\right]+0,27d", "10.3")
    m.p("em que a é a largura do guia retangular equivalente, W a distância entre os centros das fileiras e "
        "d o diâmetro das vias. Para que a fileira se comporte como parede contínua, o diâmetro e o passo s "
        "entre vias têm limites (Eqs. 6.32–6.33, s ≤ 2d). O LTA calcula a banda monomodo (de f_c a 2f_c), "
        "β, λ_g, a impedância de onda Z_TE e a perda dielétrica, e gera o layout com as vias e as transições "
        "em cunha para microstrip. A largura W concorda com a fórmula de Deslandes e Wu com diferença de "
        "0,3 a 0,6%.")
    m.figura(IMG / "tela_35_siw_estrutura.png", "10.3",
             "Aba SIW em RO3010 com corte do TE₁₀ em 6 GHz: guia com fileiras de vias e transições em cunha "
             "para microstrip de 50 Ω.")
    m.p("As transições são dimensões iniciais, a otimizar em EM. Os filtros e as linhas CRLH em SIW (seções "
        "6.5.1–6.5.2) dependem de simulação EM e ficam fora do escopo.")

    m.h("10.4  Linhas com capacitores negativos (não Foster)", 2)
    m.p("Capacitores negativos, sintetizados com circuitos ativos (conversores de impedância negativa), "
        "violam o teorema de reatância de Foster. Carregando uma linha com eles em derivação, a capacitância "
        "efetiva por unidade de comprimento diminui (" + ref("6.4.2") + "):")
    m.eq(r"C'_{\mathrm{ef}}=C'-\frac{C_{n}}{d}>0,\qquad Z=\sqrt{\frac{L'}{C'_{\mathrm{ef}}}},\qquad"
         r"\varepsilon=c^{2}L'C'_{\mathrm{ef}}<1", "10.4")
    m.p("A linha passa a ter permissividade equivalente próxima de zero (ENZ) e velocidades de fase e de "
        "grupo maiores que c, sem a dispersão das linhas passivas. A condição C'_ef > 0 é a de estabilidade. "
        "Com os dados da Fig. 6.32, o modelo dá 48,3 Ω e ε = 0,43, enquanto o livro cita 50 Ω e 0,3.")
    m.figura(IMG / "tela_36_naofoster.png", "10.4",
             "Aba Linha com capacitores negativos (Fig. 6.32): velocidades de fase e de grupo de cerca de "
             "1,53c e parâmetros S de 10 células com portas de 48,3 Ω.")
