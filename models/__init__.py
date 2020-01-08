import sys
import os
from . import ifs_model_read

thispath = os.path.dirname(os.path.abspath(__file__))

# load models
qsarnames = ['fhlb',
             'hhlb',
             'hhlt',
             'dsm',
             'tm',
             'E',
             'S',
             'A',
             'B',
             'L',
             'V',
             ]
if hasattr(sys, '_MEIPASS'):
    qsarmodels = [ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_fhlb_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlb_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlt_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_dsm_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_tm_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__E_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__S_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__A_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__B_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__L_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(sys._MEIPASS, 'ifs_qsar_V.txt')),
                  ]
else:
    qsarmodels = [ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_fhlb_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_hhlb_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_hhlt_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_dsm_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_tm_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__E_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__S_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__A_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__B_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_ADB_UFZ__L_linr.txt')),
                  ifs_model_read.QSARModel(os.path.join(thispath, 'ifs_qsar_V.txt')),
                  ]
