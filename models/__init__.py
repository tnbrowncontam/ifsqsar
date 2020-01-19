import sys
import os
from . import ifs_model_read

thispath = os.path.dirname(os.path.abspath(__file__))

if hasattr(sys, '_MEIPASS'):
    fhlb = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_fhlb_linr.txt'), 'fhlb')
    hhlb = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlb_linr.txt'), 'hhlb')
    hhlt = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlt_linr.txt'), 'hhlt')
    dsm = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_dsm_linr.txt'), 'dsm')
    tm = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_tm_linr.txt'), 'tm')
    E = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__E_linr.txt'), 'E')
    S = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__S_linr.txt'), 'S')
    A = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__A_linr.txt'), 'A')
    B = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__B_linr.txt'), 'B')
    L = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__L_linr.txt'), 'L')
    V = ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_V.txt'), 'V')
else:
    fhlb = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_fhlb_linr.txt'), 'fhlb')
    hhlb = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_hhlb_linr.txt'), 'hhlb')
    hhlt = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_hhlt_linr.txt'), 'hhlt')
    dsm = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_dsm_linr.txt'), 'dsm')
    tm = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_tm_linr.txt'), 'tm')
    E = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__E_linr.txt'), 'E')
    S = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__S_linr.txt'), 'S')
    A = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__A_linr.txt'), 'A')
    B = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__B_linr.txt'), 'B')
    L = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__L_linr.txt'), 'L')
    V = ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_V.txt'), 'V')

qsarnamelist = ['fhlb', 'hhlb', 'hhlt', 'dsm', 'tm', 'E', 'S', 'A', 'B', 'L', 'V']
qsarlist = [fhlb, hhlb, hhlt, dsm, tm, E, S, A, B, L, V]
