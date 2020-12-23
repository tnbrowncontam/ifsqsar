from . import ifs_model_read

fhlb = ifs_model_read.QSARModel('ifs_qsar_fhlb_linr', 'fhlb')
hhlb = ifs_model_read.QSARModel('ifs_qsar_hhlb_linr', 'hhlb')
hhlt = ifs_model_read.QSARModel('ifs_qsar_hhlt_linr', 'hhlt')
dsm = ifs_model_read.QSARModel('ifs_qsar_dsm_linr', 'dsm')
tm = ifs_model_read.QSARModel('ifs_qsar_tm_linr', 'tm')
E = ifs_model_read.QSARModel('ifs_qsar_ADB_UFZ__E_linr', 'E')
S = ifs_model_read.QSARModel('ifs_qsar_ADB_UFZ__S_linr', 'S')
A = ifs_model_read.QSARModel('ifs_qsar_ADB_UFZ__A_linr', 'A')
B = ifs_model_read.QSARModel('ifs_qsar_ADB_UFZ__B_linr', 'B')
L = ifs_model_read.QSARModel('ifs_qsar_ADB_UFZ__L_linr', 'L')
V = ifs_model_read.QSARModel('ifs_qsar_V', 'V')

qsarnamelist = ['fhlb', 'hhlb', 'hhlt', 'dsm', 'tm', 'E', 'S', 'A', 'B', 'L', 'V']
qsarlist = [fhlb, hhlb, hhlt, dsm, tm, E, S, A, B, L, V]
