from . import ifs_model_read

# instantiate qsar models
fhlb = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_fhlb_linr', 'fhlb')
hhlb = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_hhlb_linr', 'hhlb')
hhlt = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_hhlt_linr', 'hhlt')
dsm = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_dsm_linr', 'dsm')
tm = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_tm_linr', 'tm')
E = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__E_linr', 'E')
S = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__S_linr', 'S')
A = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__A_linr', 'A')
B = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__B_linr', 'B')
L = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_ADB_UFZ__L_linr', 'L')
V = ifs_model_read.QSARModel('ifsqsar.models.ifs_qsar_V', 'V')
logKow = ifs_model_read.METAQSARModel('ifsqsar.models.meta_qsar_logkow_pplfer', 'logKow')

# compile lists of qsar models
qsarnamelist = ['fhlb', 'hhlb', 'hhlt', 'dsm', 'tm', 'E', 'S', 'A', 'B', 'L', 'V', 'logKow']
qsarlist = [fhlb, hhlb, hhlt, dsm, tm, E, S, A, B, L, V, logKow]
