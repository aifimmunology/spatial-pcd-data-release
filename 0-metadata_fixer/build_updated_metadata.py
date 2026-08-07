"""
Build updated_metadata.csv:
  - One row per h5 file
  - Metadata joined from cohort_metadata_all.xlsx (Oncology sheet)
  - Flag files whose TIS ID is absent from the cohort sheet, or whose tissue type is non-PCD
"""
import re, os
import pandas as pd

# ── file list ──────────────────────────────────────────────────────────────
H5_FILES = [
    '../../input/3974417942/PCD_Spatial/02b58783-6669-4199-926f-e96b2b1865c6/AIFI-2026-07-01T00:14:35.136672534Z/EXP01109_Multiplex1_5_TIS05684_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/03299c70-39c4-430a-b711-dec7591eda6e/AIFI-2026-06-30T23:23:49.56310473Z/Exp998_8_chronic_myleogenous_leukemia_TIS05398-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/03ac52d2-a76a-4419-9481-b60e7b44ccbe/AIFI-2026-07-01T04:41:26.608351402Z/EXP01387_3_TIS09206-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/052a69d7-f485-4319-864d-f9c4b3a4a432/AIFI-2026-07-01T00:10:54.028781343Z/EXP01109_Multiplex1_1_TIS06392_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/0b3279f7-25cd-4655-ad9b-8ebc4879972a/AIFI-2026-06-30T23:58:32.891556606Z/EXP01109_TIS05681-001-004_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/0d2593c2-5602-48c4-bbb4-b8832f3f0b69/AIFI-2026-07-01T00:37:55.042986594Z/EXP00942_New_CD_AR_TIS06615-001-002_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/0d58fd4b-1dbd-49a9-a8a1-de6172d62e4a/AIFI-2026-06-30T23:22:43.854765602Z/Exp980_No_Stain1_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/116de923-0ca2-4042-b503-57a03742be37/AIFI-2026-07-01T00:35:37.466575988Z/EXP00904_w4_CD_AR_OldBlock_TIS06612-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/163d3b46-82d3-4d67-b4d2-d4510e20a6ed/AIFI-2026-06-30T23:23:02.612269497Z/Exp980_VIS_Stain1_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/165777a7-6d81-4333-a63d-245922d06e7e/AIFI-2026-07-01T00:35:28.242278709Z/EXP00904_w3_Visium_OldBlock_TIS06612-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/1b0bb7ef-a025-4062-86b0-bc018e318de1/AIFI-2026-07-01T04:42:53.749008921Z/EXP01387_7_TIS09210-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/1b89a1d1-e5e4-4f9c-be47-ec5cc61951f5/AIFI-2026-07-01T00:18:51.338284343Z/EXP01109_Multiplex2_5_TIS05684_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/1c48a8ba-49da-41e5-bd09-3b6940681e2c/AIFI-2026-06-30T23:58:17.386843336Z/EXP01109_TIS05681-001-004_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/1cabe22e-dc2c-4210-99b3-e412662d7f17/AIFI-2026-07-01T00:13:32.531080277Z/EXP01109_Multiplex1_4_TIS05683_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/1de9f7d6-dc27-440a-b46e-070168bce255/AIFI-2026-07-01T00:19:28.309269109Z/EXP01109_Multiplex2_5_TIS05684_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/2065ef2c-7733-46f8-81ff-589d51813a12/AIFI-2026-06-30T23:20:18.340544418Z/Exp00793_w4_No_AR_Vis_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/2143f0ad-0d91-494b-8631-598872506b20/AIFI-2026-07-01T00:38:02.886711444Z/EXP00942_Old_CD_AR_TIS06610-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/26e4e2a7-3cd8-478a-85c1-60c8a222a7b0/AIFI-2026-06-30T23:59:37.05639474Z/EXP01109_TIS05684-001-004_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/2d518247-e61d-49bd-8ddf-32ac6d0b0aba/AIFI-2026-07-01T00:13:01.266716525Z/EXP01109_Multiplex1_3_TIS05681_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/2ff8c12d-21d3-40a0-8b99-8812a80b5980/AIFI-2026-07-01T00:31:15.943595914Z/EXP00998_1_MM_TIS05392-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/339b651d-4a55-4212-8b81-075277bffed2/AIFI-2026-07-01T00:11:26.986254444Z/EXP01109_Multiplex1_1_TIS06392_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/36146546-c397-447e-8638-5368c0fd087b/AIFI-2026-06-30T23:23:33.895564825Z/Exp998_6_early_myoproliferative_TIS05395-001-009_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/37654a93-1335-477a-af3d-372cd3303d9b/AIFI-2026-06-30T23:20:59.213736463Z/Exp00904_w2_CD_AR_NewBlock_TIS05740-001-018_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/3e603440-2dab-4963-a191-b231ee5e08a4/AIFI-2026-06-30T23:23:18.302983041Z/Exp998_10_MM_Plasmacytoma_TIS06392-001-003_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/3ed83e9b-c145-403f-8311-4b75f88ec004/AIFI-2026-07-01T04:41:03.19250372Z/EXP01387_2_TIS09205-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/46cfe85a-ac73-408d-8c3d-b4f234251113/AIFI-2026-06-30T23:22:20.419002879Z/Exp942_New_Vis_AR_TIS06615-001-002_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/486dfc1d-003f-4c67-8477-61187a79e703/AIFI-2026-06-30T23:23:41.817197928Z/Exp998_7_MM_TIS05395-001-009_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/4ad23f67-7fa7-4fbe-a0fe-a54b8ab14145/AIFI-2026-07-01T00:31:23.943820519Z/EXP00998_3_MM_TIS05397-001-009_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/4bfa07c7-9b52-4821-903d-e644f2269f7f/AIFI-2026-07-01T00:01:14.143750482Z/EXP01109_TIS06392-002-005_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/4e740861-1d06-478e-8390-f6ec5774714e/AIFI-2026-07-01T00:15:06.140527866Z/EXP01109_Multiplex2_1_TIS06392_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/4f6dfefe-cad7-46a4-9074-21576d32e599/AIFI-2026-07-01T00:15:37.384443888Z/EXP01109_Multiplex2_2_TIS06377_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/516bb82e-328c-4ec6-a054-e3ea8855494b/AIFI-2026-06-30T23:22:36.032810495Z/Exp980_CST_Stain1_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/567d5b34-3ec0-4c64-b986-f760a09c55fb/AIFI-2026-06-30T23:59:21.433285898Z/EXP01109_TIS05684-001-004_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/5e2005d8-6e8c-4ce5-ba94-50cf01b4c9cc/AIFI-2026-07-01T00:31:39.650364179Z/EXP00998_5_MM_TIS05393-001-010_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/5ef32b92-e360-429d-9154-d496da66935f/AIFI-2026-06-30T23:20:34.177114609Z/Exp00793_w6_AR_No_Staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/5fe7a47b-11c4-438f-a262-71cf86cae152/AIFI-2026-07-01T04:42:30.326545791Z/EXP01387_6_TIS09209-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/6267c044-e9ee-4297-a1a9-65fddecb257c/AIFI-2026-06-30T23:20:10.442772498Z/Exp00793_w3_No_AR_Vis_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/62bb3015-d015-4eeb-b17e-390ef82473c3/AIFI-2026-07-01T04:42:05.609579046Z/EXP01387_5_TIS09208-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/64271569-8bd8-4632-8dc2-123fdbef6bf0/AIFI-2026-07-01T00:00:58.456808541Z/EXP01109_TIS06392-002-005_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/66a68627-0daf-4021-a565-60d76559a0f5/AIFI-2026-07-01T00:35:18.650627389Z/EXP00904_w1_Visium_NewBlock_TIS05740-001-018_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/693df872-150a-4f9e-83bc-ab448cc4b987/AIFI-2026-07-01T00:12:14.059380758Z/EXP01109_Multiplex1_2_TIS06377_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/6b6bd4d0-6a38-42bf-9f3c-f56d493f8af3/AIFI-2026-07-01T00:12:45.503916719Z/EXP01109_Multiplex1_3_TIS05681_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/6c41e6bb-6fb4-4636-8a84-532c2912eae1/AIFI-2026-06-30T23:22:04.882829763Z/Exp01349_TissDiss_SP7_TIS06384-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/705fc4ac-0dd1-4ad8-991f-b98df1be0ad3/AIFI-2026-07-01T04:43:09.571367017Z/EXP01387_8_TIS09211-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/7089ed8e-45fc-444b-afb8-ec6ee631cc75/AIFI-2026-07-01T00:40:02.261468932Z/EXP00980_Dub_Stain2_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/7526ebf9-aa9c-4f3b-ba58-87535adb34f7/AIFI-2026-06-30T23:21:33.673641731Z/Exp01349_TissDiss_SP3_TIS06665-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/78deb36e-bd93-4f39-b425-f62f3c90db58/AIFI-2026-07-01T00:13:16.897191822Z/EXP01109_Multiplex1_4_TIS05683_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/7986e584-c1d7-4a9e-86e3-34037681430b/AIFI-2026-07-01T00:01:45.498975665Z/EXP01109_TIS06402-001-005_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/79cfb7f1-527c-457e-9052-c56b679dbaa4/AIFI-2026-06-30T23:22:12.573191812Z/Exp01349_TissDiss_SP8_TIS06397-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/7e26ca5f-0bfa-418e-a93f-8e3de7de6cc0/AIFI-2026-06-30T23:22:28.26044787Z/Exp942_Old_Vis_AR_TIS06610-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/81af423d-6480-4b73-8bd6-c1a8eb54abe9/AIFI-2026-07-01T00:13:48.284552123Z/EXP01109_Multiplex1_4_TIS05683_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/83a08a14-f58a-49e6-9ade-7cc76d79d1fc/AIFI-2026-06-30T23:21:49.446010622Z/Exp01349_TissDiss_SP5_TIS06394-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/86377f12-3fb1-4777-84b2-c4d213347fb9/AIFI-2026-07-01T00:15:21.801130721Z/EXP01109_Multiplex2_1_TIS06392_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/895e7fb4-0949-489b-9786-9bea6f550138/AIFI-2026-07-01T00:17:52.189967135Z/EXP01109_Multiplex2_4_TIS05683_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/8b5ba6f2-9784-4402-9e7f-96d8da3b7179/AIFI-2026-06-30T23:59:52.745193422Z/EXP01109_TIS06377-001-005_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/91970e83-102b-48a2-b772-0b90cf4797ef/AIFI-2026-07-01T00:11:42.470044333Z/EXP01109_Multiplex1_2_TIS06377_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/93420983-831a-4060-979d-5385ede93a02/AIFI-2026-07-01T00:14:19.54765505Z/EXP01109_Multiplex1_5_TIS05684_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/93db867c-e5ef-4bb6-ac0d-3f2d69d4515e/AIFI-2026-07-01T00:14:04.043195753Z/EXP01109_Multiplex1_5_TIS05684_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/953d4213-599b-422e-81dc-8f6dd65653b6/AIFI-2026-07-01T00:31:31.762745093Z/EXP00998_4_Normal_TIS04788-001-011_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/9579044a-3bd2-413a-9e76-765f26d12bff/AIFI-2026-07-01T00:17:24.850233351Z/EXP01109_Multiplex2_4_TIS05683_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/9a0d139c-0859-4ae0-beb5-52f21665173c/AIFI-2026-07-01T00:12:29.927731396Z/EXP01109_Multiplex1_3_TIS05681_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/9b0accf5-0f6d-475f-b6ac-809ed171d873/AIFI-2026-07-01T00:16:40.044247521Z/EXP01109_Multiplex2_3_TIS05681_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/a5101e80-656a-4ac2-9bfc-f955b4879b2d/AIFI-2026-07-01T00:14:50.691416705Z/EXP01109_Multiplex2_1_TIS06392_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/a55560fb-47a5-416e-aca4-7d4570aefabc/AIFI-2026-06-30T23:19:54.542611682Z/Exp00793_w1_No_AR_No_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/aace7f2f-a03c-45be-a43b-6043624e0565/AIFI-2026-07-01T04:40:39.342561158Z/EXP01387_1_TIS09204-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/ab36a3e8-e025-4088-a65b-48413d14b840/AIFI-2026-07-01T00:15:53.095081523Z/EXP01109_Multiplex2_2_TIS06377_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/ae5fe019-dd17-4d22-91c0-205a3cb26a4a/AIFI-2026-06-30T23:59:04.252941227Z/EXP01109_TIS05683-001-004_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/b0981942-ee8b-47bd-bc7e-c5c8ce4b7ccd/AIFI-2026-07-01T00:01:29.894522876Z/EXP01109_TIS06402-001-005_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/b33aa5cf-6d85-40ca-9206-4323e54ad05b/AIFI-2026-07-01T04:41:42.130363037Z/EXP01387_4_TIS09207-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/b67736dc-528d-4d36-9b71-c5c91ded75cc/AIFI-2026-06-30T23:20:26.258017378Z/Exp00793_w5_AR_No_Staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/b744c25c-6d2a-4ce6-be0e-8ef90cab928c/AIFI-2026-06-30T23:21:41.433558541Z/Exp01349_TissDiss_SP4_TIS06664-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/bb78b058-2dc5-42f9-b5c7-e196ddd5e4f1/AIFI-2026-06-30T23:21:18.00831889Z/Exp01349_TissDiss_SP1_TIS06385-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/beac980f-e475-4e95-ac2a-f1766c114db2/AIFI-2026-07-01T00:11:58.048712399Z/EXP01109_Multiplex1_2_TIS06377_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c02a952c-612a-43fa-a251-67dbf2f83bda/AIFI-2026-07-01T00:16:56.910076171Z/EXP01109_Multiplex2_3_TIS05681_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c1d76fbc-ba00-432b-9b7a-48f91cffd0bc/AIFI-2026-06-30T23:20:02.685444474Z/Exp00793_w2_No_AR_No_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c2898b36-3bed-44e9-a28c-708033829548/AIFI-2026-07-01T00:00:10.568789255Z/EXP01109_TIS06377-001-005_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c59067e6-2973-4151-90ad-529fdd804403/AIFI-2026-06-30T23:21:25.784012568Z/Exp01349_TissDiss_SP2_TIS06389-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c6c993d6-c758-458c-b918-1e4fe11db0ba/AIFI-2026-07-01T00:16:24.37188031Z/EXP01109_Multiplex2_3_TIS05681_1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c71f2420-9a1b-4486-ae43-cae08e737464/AIFI-2026-07-01T00:39:54.502926606Z/EXP00980_CST_Stain2_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/c86be991-caee-4976-a528-29144501b267/AIFI-2026-07-01T00:16:08.806824273Z/EXP01109_Multiplex2_2_TIS06377_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/d233790d-8242-4ece-be3e-f8b37ee7b2fd/AIFI-2026-07-01T00:11:09.735916274Z/EXP01109_Multiplex1_1_TIS06392_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/d565e68d-0a0b-4661-9fa6-25798b6d855f/AIFI-2026-07-01T00:31:47.745789536Z/EXP00998_9_MM_TIS04789-001-015_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/d6c3b7a9-1398-44b1-8abc-9045c8fe34f6/AIFI-2026-07-01T00:00:27.273899623Z/EXP01109_TIS06381-002-005_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/e2835f08-8510-4b76-bcd1-575fe44ba6c3/AIFI-2026-06-30T23:23:26.049519835Z/Exp998_2_Normal_TIS05396-001-009_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/e58ae0cb-e635-4d7e-a311-32ec1bc32cca/AIFI-2026-06-30T23:58:48.4152182Z/EXP01109_TIS05683-001-004_c1_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/ee0457b6-7cbf-40f8-9d33-66bfeec8405a/AIFI-2026-07-01T00:00:42.974880058Z/EXP01109_TIS06381-002-005_c2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/ef3fee0a-4b3a-4280-bdf7-a18f3a4ca2a8/AIFI-2026-06-30T23:22:51.867713414Z/Exp980_No_Stain2_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/efd370be-7277-4992-82bf-e9c19196b668/AIFI-2026-07-01T00:19:12.717422386Z/EXP01109_Multiplex2_5_TIS05684_2_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/f04e6cfe-2567-4c2c-9e5b-1ccdd7e4a600/AIFI-2026-06-30T23:20:49.678512911Z/Exp00793_w8_AR_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/f131e63e-e536-4399-8a04-6c5d22a18371/AIFI-2026-07-01T00:18:20.75724557Z/EXP01109_Multiplex2_4_TIS05683_3_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/f2e82098-1980-4fc6-a7fd-213caf4422ee/AIFI-2026-06-30T23:23:10.558264254Z/Exp980_VIS_Stain2_TIS05841-001-008_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/fa87eb82-61e5-4e6f-8b75-d9f7c39a114a/AIFI-2026-06-30T23:20:41.884757517Z/Exp00793_w7_AR_staining_TIS05740-001-001_sample_filtered_feature_bc_matrix.h5',
    '../../input/3974417942/PCD_Spatial/fc91fc1f-db57-4887-a6ae-eef70c52f537/AIFI-2026-06-30T23:21:57.117695089Z/Exp01349_TissDiss_SP6_TIS06390-001_sample_filtered_feature_bc_matrix.h5',
]

# ── helpers (same logic as notebook) ───────────────────────────────────────
def extract_tis(fname):
    m = re.search(r'(TIS\d{5}(?:-\d{3}(?:-\d{3})?)?)', fname)
    return m.group(1) if m else 'unknown'

def infer_project(fname):
    f = fname.upper()
    if 'EXP00998' in f or 'EXP998' in f:       return 'aacr1'
    if 'EXP01109' in f and 'MULTIPLEX' in f:    return 'aacr2-mp'
    if 'EXP01109' in f:                         return 'aacr2-sp'
    if 'EXP01349' in f and 'MULTIPLEX' in f:    return 'aacr3-mp'
    if 'EXP01349' in f:                         return 'aacr3-sp'
    if 'EXP01387' in f:                         return 'aacr4'
    return 'unknown'

def tis_base(tis_id):
    """Return the Oncology-sheet lookup key for a TIS ID.
    TIS05684-001-004 → TIS05684-001
    TIS05684-001     → TIS05684-001
    TIS05684         → TIS05684-001  (Multiplex filenames omit the -001 suffix)
    """
    m = re.match(r'(TIS\d{5}-\d{3})', tis_id)
    if m:
        return m.group(1)
    # bare 5-digit code: assume first accession block (-001)
    m5 = re.match(r'TIS\d{5}', tis_id)
    return m5.group() + '-001' if m5 else tis_id

# ── build file-level table ─────────────────────────────────────────────────
rows = []
for path in H5_FILES:
    fname   = os.path.basename(path)
    src     = fname.replace('_sample_filtered_feature_bc_matrix.h5', '')
    tis_id  = extract_tis(fname)
    project = infer_project(fname)
    rows.append({'file_path': path, 'source_file': src,
                 'TIS_ID': tis_id, 'project': project})

df = pd.DataFrame(rows)
df['TIS_ID_base'] = df['TIS_ID'].apply(tis_base)

# ── load cohort metadata ───────────────────────────────────────────────────
oc = pd.read_excel('cohort_metadata_all.xlsx', sheet_name='Oncology')
keep = ['AIFI Barcode', 'Subject ID', 'For Experiments', 'Tissue Type',
        'Sample Site', 'Disease State', 'Patient Age at collection', 'Sex',
        'Race/Ethnicity', 'Treatment (if known)', '%Tumor',
        'Collection Year', 'Biopsy or Resection']
oc_clean = (oc[keep]
            .dropna(subset=['AIFI Barcode'])
            .rename(columns={
                'AIFI Barcode':              'TIS_ID_base',
                'Subject ID':               'subject_id',
                'For Experiments':          'cohort_label',
                'Tissue Type':              'tissue_type',
                'Sample Site':              'sample_site',
                'Disease State':            'disease_state',
                'Patient Age at collection':'age_at_collection',
                'Sex':                      'sex',
                'Race/Ethnicity':           'race_ethnicity',
                'Treatment (if known)':     'treatment',
                '%Tumor':                   'pct_tumor',
                'Collection Year':          'collection_year',
                'Biopsy or Resection':      'biopsy_or_resection',
            }))

# ── join ───────────────────────────────────────────────────────────────────
merged = df.merge(oc_clean, on='TIS_ID_base', how='left')

# ── flags ──────────────────────────────────────────────────────────────────
NON_PCD_TISSUES = {'tonsil', 'lymph node'}

def make_flag(row):
    if pd.isna(row['tissue_type']) and pd.isna(row['disease_state']) and pd.isna(row['subject_id']):
        return 'NOT_IN_COHORT_METADATA'
    tt = str(row['tissue_type']).strip().lower() if not pd.isna(row['tissue_type']) else ''
    if tt in NON_PCD_TISSUES:
        return f'NON_PCD_TISSUE:{row["tissue_type"]}'
    ds = str(row['disease_state']).strip().lower() if not pd.isna(row['disease_state']) else ''
    if ds in ('normal', 'nan', ''):
        # check explicitly 'Normal BM' label in cohort_label
        cl = str(row.get('cohort_label', '')).lower()
        if 'normal' in cl:
            return 'NORMAL_CONTROL'
    return ''

merged['flag'] = merged.apply(make_flag, axis=1)
merged['in_cohort_metadata'] = merged['flag'] != 'NOT_IN_COHORT_METADATA'

# reorder: identity cols first, flags, then metadata
col_order = [
    'file_path', 'source_file', 'TIS_ID', 'TIS_ID_base', 'project',
    'in_cohort_metadata', 'flag',
    'subject_id', 'cohort_label', 'tissue_type', 'sample_site',
    'disease_state', 'age_at_collection', 'sex', 'race_ethnicity',
    'treatment', 'pct_tumor', 'collection_year', 'biopsy_or_resection',
]
merged = merged[col_order]
merged.to_csv('updated_metadata.csv', index=False)
print(merged[['source_file','TIS_ID_base','project','in_cohort_metadata','flag']].to_string())
print(f"\n{len(merged)} files | flags: {merged['flag'].value_counts().to_dict()}")
