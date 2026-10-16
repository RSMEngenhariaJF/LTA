# LTA – Linhas de Transmissão Artificiais Periódicas

Ferramenta de projeto com interface gráfica para linhas periódicas (F. Martín, 2015, cap. 2): célula ABCD, linha com carga capacitiva (slow-wave), EBG por modos acoplados, análise de Floquet (harmônicos espaciais e leaky-wave), linha de atraso chirpada e ressonador de Bragg, com perdas e dispersão da microstrip opcionais. Do cap. 3: síntese de linhas CRLH/D-CRLH, linha microstrip CRLH com interdigitais e stubs com via, linha com CSRR e gap, e extração de parâmetros de circuito a partir de .s2p (Apêndice G). Exporta para KiCad, Altium (DXF/Gerber/KiCad), ADS (GDSII/DXF/Touchstone), HFSS (script) e CST (macro).

```
pip install -r requirements.txt
python main.py                          # ou executar.bat
python -m unittest discover testes -v   # validação contra o livro
```

- Manual completo: [doc/01_manual.md](doc/01_manual.md)
- Roteiros: [cap. 2](doc/00_roteiro_cap2.md), [cap. 3](doc/00_roteiro_cap3.md)

Os resultados são estimativas de modelos de circuito e quase-estáticos. Valide por simulação EM e medição, com revisão técnica, antes de usar.
