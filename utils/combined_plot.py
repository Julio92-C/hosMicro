import PyPDF2


merger = PyPDF2.PdfFileMerger()

file_names = ["./Species_count_abricate_Kraken2.pdf",
              "./AMR_count_abricate_Kraken2.pdf",
              "./VFs_count_abricate_Kraken2.pdf",
              "./Plasmid_count_abricate_Kraken2.pdf",
              "./AMR_VFs_MGEs_abri_Kraken2_stacked_barplot.pdf",
              "./AMR_VFs_MGEs_abri_Kraken2_pieChart.pdf",
              "./core_accesoryAMR_VFs_MGEs_abricate_Kraken2_GeomBin.pdf",
              "./coreAMR_VFs_MGEs_abri_Kraken2Update1-95.pdf",
              "Pathogen_AMR_VFs_MGEs_all_pheatmap.pdf",
              "Pathogen_AMR_VFs_MGEs_100rows_pheatmap.pdf",
              "Air_Sur_Venndiagram.pdf",
              "AMR_VFs_MGEs_vennDiagram.pdf"
              ]


for file_name in file_names:
    merger.append(file_name)

merger.write("./Species_AMR_VFs_MGEs_basic_stats.pdf")
merger.close()
