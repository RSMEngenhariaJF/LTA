# LTA – Linhas de Transmissão Artificiais Periódicas

Ferramenta de projeto com interface gráfica para linhas periódicas (F. Martín, 2015, cap. 2): célula ABCD, linha com carga capacitiva (slow-wave) e EBG por modos acoplados. Exporta para KiCad, Altium (DXF/Gerber/KiCad), ADS (GDSII/DXF/Touchstone), HFSS (script) e CST (macro).

```
pip install -r requirements.txt
python main.py                          # ou executar.bat
python -m unittest discover testes -v   # validação contra o livro
```

- Manual completo: [doc/01_manual.md](doc/01_manual.md)
- Roteiro do capítulo: [doc/00_roteiro_cap2.md](doc/00_roteiro_cap2.md)

Os resultados são estimativas de modelos de circuito e quase-estáticos. Valide por simulação EM e medição, com revisão técnica, antes de usar.
