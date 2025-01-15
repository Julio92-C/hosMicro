import PyPDF2


merger = PyPDF2.PdfFileMerger()

file_names = ["./taxaCount_barPlot.pdf",
              "./airSur_vennDiagram.pdf",
              "./countTaxa_violinPlot.pdf",
              "./relativeAbundance_stackedPlot1.pdf",
              "./countARG_barPlot.pdf",
              "./countARG_violinPlot.pdf",
              "./countVFs_barPlot.pdf",
              "./countVFs_violinPlot.pdf",
              "./countMGEs_barPlot.pdf",
              "./disGEs_stackedBarplot.pdf",
              "./geneCount_pieChart.pdf",
              "./geneticElements_vennDiagram.pdf",
              "./geneticElem_PCoA.pdf",
              "./airSur_chordDiagram.pdf",
              "./airMicrobiome_pheatmap.pdf",
              ]


for file_name in file_names:
    merger.append(file_name)

merger.write("./HosMicro_Summary.pdf")
merger.close()
