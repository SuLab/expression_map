library(recount)

tcga.meta <- all_metadata(subset='tcga')

tcga.meta <- tcga.meta[,c(19,1,66,78,75,79,92,110,94,93,68,71,41,87,29,198)]

colnames(tcga.meta) <- c('data_id','project','sample_id','tissue_general','tissue_detail','tumor_id','tumor_stage','tumor_sample_type','living','age','sex','ethnicity','center_submittter','center_source','seq_platform','treatment_result')

tcga.meta$organ <- NA
tcga.meta$tissue <- NA
tcga.meta$efo_disease <- NA
tcga.meta$doid_disease <- NA

## Liver

tcga.meta$organ[tcga.meta$tissue_detail == 'Liver Hepatocellular Carcinoma'] <- '<UBERON:0002107: liver>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Liver Hepatocellular Carcinoma'] <- '<UBERON:0002107: liver>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Liver Hepatocellular Carcinoma'] <- '<EFO:0000182: hepatocellular carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Liver Hepatocellular Carcinoma'] <- '<DOID:684: hepatocellular carcinoma>'

## Pleura

tcga.meta$organ[tcga.meta$tissue_general == 'Pleura'] <- '<UBERON:0000977: pleura>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Pleura'] <- '<UBERON:0003390: mesothelium of pleural cavity>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Pleura'] <- '<EFO:0000588: mesothelioma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Pleura'] <- '<DOID:1790: mesothelioma>'

## Stomach

tcga.meta$organ[tcga.meta$tissue_general == 'Stomach'] <- '<UBERON:0000945: stomach>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Stomach'] <- '<UBERON:0001276: epithelium of stomach>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Stomach'] <- '<EFO:0000503: gastric adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Stomach'] <- '<DOID:3717: gastric adenocarcinoma>'

## Adrenal Gland: paraganglioma

tcga.meta$organ[tcga.meta$tissue_general == 'Adrenal Gland'] <- '<UBERON:0002369: adrenal gland>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Pheochromocytoma and Paraganglioma'] <- '<UBERON:0002369: adrenal gland>'

######### unsure about this... need to decide on paraganlioma / pheochromocytoma OR just adrenal cancer

tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Pheochromocytoma and Paraganglioma'] <- '<EFO:0000182: hepatocellular carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Pheochromocytoma and Paraganglioma'] <- '<DOID:684: hepatocellular carcinoma>'

## Adrenal Gland: adrenocortical carcinoma

tcga.meta$tissue[tcga.meta$tissue_detail == 'Adrenocortical Carcinoma'] <- '<UBERON:0001235: adrenal cortex>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Adrenocortical Carcinoma'] <- '<EFO:0003093: adrenocortical carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Adrenocortical Carcinoma'] <- '<DOID:3948: adrenocortical carcinoma>'

## Pancreas

tcga.meta$organ[tcga.meta$tissue_general == 'Pancreas'] <- '<UBERON:0001264: pancreas>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Pancreas'] <- '<UBERON:0009970: epithelium of pancreatic duct>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Pancreas'] <- '<EFO:1000044: pancreatic adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Pancreas'] <- '<DOID:4074: pancreatic adenocarcinoma>'

## Uterine Corpus Edometrial Carcinoma

tcga.meta$organ[tcga.meta$tissue_detail == 'Uterine Corpus Endometrial Carcinoma'] <- '<UBERON:0000995: uterus>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Uterine Corpus Endometrial Carcinoma'] <- '<UBERON:0001295: endometrium>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Uterine Corpus Endometrial Carcinoma'] <- '<EFO:1001512: endometrial carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Uterine Corpus Endometrial Carcinoma'] <- '<DOID:0050939: uterine corpus endometrial carcinoma>'

## Uterine carcinosarcoma

tcga.meta$organ[tcga.meta$tissue_detail == 'Uterine Carcinosarcoma'] <- '<UBERON:0000995: uterus>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Uterine Carcinosarcoma'] <- '<UBERON:0000459: uterine wall>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Uterine Carcinosarcoma'] <- '<EFO:1000613: uterine carcinosarcoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Uterine Carcinosarcoma'] <- '<DOID:6171: uterine carcinosarcoma>'

## Testis

tcga.meta$organ[tcga.meta$tissue_general == 'Testis'] <- '<UBERON:0000473: testis>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Testis'] <- '<UBERON:0000473: testis>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Testis'] <- '<EFO:1000566: testicular germ cell tumor>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Testis'] <- '<DOID:5557: testicular germ cell cancer>'

## Kidney Rencal clear Cell Carcinoma

tcga.meta$organ[tcga.meta$tissue_general == 'Kidney'] <- '<UBERON:0002113: kidney>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Kidney Renal Clear Cell Carcinoma'] <- '<UBERON:0004819: kidney epithelium>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Kidney Renal Clear Cell Carcinoma'] <- '<EFO:0000349: clear cell renal carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Kidney Renal Clear Cell Carcinoma'] <- '<DOID:4467: renal clear cell carcinoma>'

## Kidney Chromophobe

tcga.meta$tissue[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<UBERON:0001285: nephron>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<EFO:0000335: chromophobe renal cell carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<DOID:4471: chromophobe renal cell carcinoma>'

## Kidney Renal Papillary Cell Carcinoma

tcga.meta$tissue[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<UBERON:0005167: papillary duct>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<EFO:0000640: papillary renal cell carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<DOID:4465: papillary renal cell carcinoma>'

## Glioblastoma Multiforme

tcga.meta$organ[tcga.meta$tissue_general == 'Brain'] <- '<UBERON:0000955: brain>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Brain'] <- '<UBERON:0003714: neural tissue>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Glioblastoma Multiforme'] <- '<EFO:0000519: glioblastoma multiforme>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Kidney Chromophobe'] <- '<DOID:3073: brain glioblastoma multiforme>'

## Glioma

tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Glioma'] <- '<EFO:0005543: glioma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Glioma'] <- '<EFO:0005543: glioma>'

## Bladder Urothelial Carcinoma

tcga.meta$organ[tcga.meta$tissue_general == 'Bladder'] <- '<UBERON:0001255: urinary bladder>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Bladder'] <- '<UBERON:0001255: urinary bladder>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Bladder'] <- '<EFO:0000292: bladder carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Bladder'] <- '<DOID:4006: bladder urothelial carcinoma>'

## Esophagus

tcga.meta$organ[tcga.meta$tissue_general == 'Esophagus'] <- '<UBERON:0001043: esophagus>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Esophagus'] <- '<UBERON:0001976: epithelium of esophagus>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Esophagus'] <- '<EFO:0002916: esophageal carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Esophagus'] <- '<DOID:1107: esophageal carcinoma>'

## Skin

tcga.meta$organ[tcga.meta$tissue_general == 'Skin'] <- '<UBERON:0000014: zone of skin>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Skin'] <- '<UBERON:0002025: stratum basale of epidermis>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Skin'] <- '<EFO:0000389: cutaneous melanoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Skin'] <- '<DOID:1909: melanoma>'

## Cervical Squamous Cell Carcinoma

tcga.meta$organ[tcga.meta$tissue_general == 'Cervix'] <- '<UBERON:0000002: uterine cervix>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Cervix'] <- '<UBERON:00004801: cervix epithelium>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Cervical Squamous Cell Carcinoma'] <- '<EFO:0000389: cutaneous melanoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Cervical Squamous Cell Carcinoma'] <- '<DOID:3744: cervical squamous cell carcinoma>'

## Endocervical Adenocarcinoma

tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Endoccervical Adenocarcinoma'] <- '<EFO:0001416: endocervical adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Endocervical Adenocarcinoma'] <- '<DOID:0050940: endocervical adenocarcinoma>'

## Bone Marrow

tcga.meta$organ[tcga.meta$tissue_general == 'Bone Marrow'] <- '<UBERON:0002371: bone marrow>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Bone Marrow'] <- '<UBERON:0002371: bone marrow>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Bone Marrow'] <- '<EFO:0000222: acute myeloid leukemia>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Bone Marrow'] <- '<DOID:9119: acute myeloid leukemia>'

## Lung

tcga.meta$organ[tcga.meta$tissue_general == 'Lung'] <- '<UBERON:0002048: lung>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Lung'] <- '<UBERON:0000115: lung epithelium>'

## Lung Adenocarcinoma

tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Lung Adenocarcinoma'] <- '<EFO:0000571: lung adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Lung Adenocarcinoma'] <- '<DOID:3910: lung adenocarcinoma>'

## Lung Squamous Cell Carcinoma

tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Lung Squamous Cell Carcinoma'] <- '<EFO:0000708: squamous cell lung carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Lung Squamous Cell Carcinoma'] <- '<DOID:3907: lung sqamous cell carcinoma>'

## Breast

tcga.meta$organ[tcga.meta$tissue_general == 'Breast'] <- '<UBERON:0000310: breast>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Breast'] <- '<UBERON:0008367: breast epithelium>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Breast'] <- '<EFO:0000305: breast carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Breast'] <- '<DOID:3459: breast carcinoma>'

## Prostate

tcga.meta$organ[tcga.meta$tissue_general == 'Prostate'] <- '<UBERON:0002367: prostate gland>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Prostate'] <- '<UBERON:0000428: prostate epithelium>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Prostate'] <- '<EFO:0000673: prostate adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Prostate'] <- '<DOID:2526: prostate adenocarcinoma>'

## Eye

tcga.meta$organ[tcga.meta$tissue_general == 'Eye'] <- '<UBERON:0000970: eye>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Eye'] <- '<UBERON:0001768: uvea>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Eye'] <- '<EFO:1000616: uveal melanoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Eye'] <- '<DOID:6039: uveal melanoma>'

## Thymus

tcga.meta$organ[tcga.meta$tissue_general == 'Thymus'] <- '<UBERON:0002370: thymus>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Thymus'] <- '<UBERON:0003846: thymus epithelium>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Thymus'] <- '<DOID:3275: thymoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Thymus'] <- '<EFO:1000581: Thymoma>'

## Thyroid

tcga.meta$organ[tcga.meta$tissue_general == 'Thyroid'] <- '<UBERON:0002046: thyroid gland>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Thyroid'] <- '<UBERON:0002046: thyroid gland>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Thyroid'] <- '<EFO:0002892: thyroid carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Thyroid'] <- '<DOID:3963: thyroid carcinoma>'

## Rectum Adenocarincoma

tcga.meta$organ[tcga.meta$tissue_detail == 'Rectum Adenocarcinoma'] <- '<UBERON:0012652: colorectum>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Rectum Adenocarcinoma'] <- '<UBERON:0001052: rectum>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Rectum Adenocarcinoma'] <- '<EFO:0005631: rectal adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Rectum Adenocarcinoma'] <- '<DOID:1996: rectal adenocarcinoma>'

## Colon Adenocarcinoma

tcga.meta$organ[tcga.meta$tissue_detail == 'Colon Adenocarcinoma'] <- '<UBERON:0012652: colorectum>'
tcga.meta$tissue[tcga.meta$tissue_detail == 'Colon Adenocarcinoma'] <- '<UBERON:0001155: colon>'
tcga.meta$efo_disease[tcga.meta$tissue_detail == 'Colon Adenocarcinoma'] <- '<EFO:1001949: colon adenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_detail == 'Colon Adenocarcinoma'] <- '<DOID:234: colon adenocarcinoma>'

## Soft Tissue

tcga.meta$organ[tcga.meta$tissue_general == 'Soft Tissue'] <- '<UBERON:0003104: mesenchyme>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Soft Tissue'] <- '<UBERON:0003104: mesenchyme>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Soft Tissue'] <- '<EFO:1001968: soft tissue sarcoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Soft Tissue'] <- '<DOID:1115: sarcoma>'

## Lymph Nodes

tcga.meta$organ[tcga.meta$tissue_general == 'Lymph Nodes'] <- '<UBERON:0000029: lymph node>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Lymph Nodes'] <- '<CL:0000236: B cell>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Lymph Nodes'] <- '<EFO:0000403: diffuse large B-cell lymphoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Lymph Nodes'] <- '<DOID:0050745: diffuse large B-cell lymphoma>'

## Ovary

tcga.meta$organ[tcga.meta$tissue_general == 'Ovary'] <- '<UBERON:0000992: female gonad>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Ovary'] <- '<CL:2000063: ovarian fibroblast>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Ovary'] <- '<EFO:1000043: ovarian serous cystadenocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Ovary'] <- '<DOID:5746: ovarian serous cystadenocarcinoma>'

## Head and Neck

tcga.meta$organ[tcga.meta$tissue_general == 'Head and Neck'] <- '<UBERON:0007811: craniocervical region>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Head and Neck'] <- '<UBERON:0007811: craniocervical region>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Head and Neck'] <- '<EFO:0000181: head and neck squamous cell carcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Head and Neck'] <- '<DOID:5520: head and neck squamous cell carcinoma>'

## Bile Duct

tcga.meta$organ[tcga.meta$tissue_general == 'Bile Duct'] <- '<UBERON:0002394: bile duct>'
tcga.meta$tissue[tcga.meta$tissue_general == 'Bile Duct'] <- '<UBERON:0002394: bile duct>'
tcga.meta$efo_disease[tcga.meta$tissue_general == 'Bile Duct'] <- '<EFO:0005221: cholangiocarcinoma>'
tcga.meta$doid_disease[tcga.meta$tissue_general == 'Bile Duct'] <- '<DOID:4947: cholangiocarcinoma>'
