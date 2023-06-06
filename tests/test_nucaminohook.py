import unittest

from sierralocal.nucaminohook import NucAminoAligner
from sierralocal.hivdb import HIVdb

class TestNucAminoAligner(unittest.TestCase):
    def setUp(self):
        self.algorithm = HIVdb()
        self.aligner = NucAminoAligner(self.algorithm, program='nuc')
        self.postAligner = NucAminoAligner(self.algorithm)

    def testPrevalenceParser(self):
        # Setting params
        filename = 'PIPrevalences.tsv'

        expected = {'1PLA': 0.0, '1PLB': 0.1, '1PLC': 0.0, '1PLD': 0.0, '1PLF': 0.0, '1PLG': 0.0, '1PLCRF01_AE': 0.0, '1PLCRF02_AG': 0.0, '1PLOther': 0.1, '1PLAll': 0.1, '1PSA': 0.0, '1PSB': 0.0, '1PSC': 0.1, '1PSD': 0.0, '1PSF': 0.0, '1PSG': 0.3, '1PSCRF01_AE': 0.0, '1PSCRF02_AG': 0.0, '1PSOther': 0.0, '1PSAll': 0.0, '1PRA': 0.4, '1PRB': 0.0, '1PRC': 0.0, '1PRD': 0.0, '1PRF': 0.0, '1PRG': 0.0, '1PRCRF01_AE': 0.0, '1PRCRF02_AG': 0.0, '1PROther': 0.0, '1PRAll': 0.0, '1PEA': 0.0, '1PEB': 0.0, '1PEC': 0.0, '1PED': 0.0, '1PEF': 0.0, '1PEG': 0.0, '1PECRF01_AE': 0.0, '1PECRF02_AG': 0.0, '1PEOther': 0.1, '1PEAll': 0.0, '1PFA': 0.0, '1PFB': 0.0, '1PFC': 0.0, '1PFD': 0.0, '1PFF': 0.0, '1PFG': 0.0, '1PFCRF01_AE': 0.0, '1PFCRF02_AG': 0.0, '1PFOther': 0.0, '1PFAll': 0.0, '1PTA': 0.0, '1PTB': 0.0, '1PTC': 0.0, '1PTD': 0.0, '1PTF': 0.0, '1PTG': 0.0, '1PTCRF01_AE': 0.0, '1PTCRF02_AG': 0.0, '1PTOther': 0.0, '1PTAll': 0.0, '1PQA': 0.0, '1PQB': 0.0, '1PQC': 0.0, '1PQD': 0.0, '1PQF': 0.0, '1PQG': 0.0, '1PQCRF01_AE': 0.0, '1PQCRF02_AG': 0.0, '1PQOther': 0.0, '1PQAll': 0.0, '1PAA': 0.0, '1PAB': 0.0, '1PAC': 0.0, '1PAD': 0.0, '1PAF': 0.0, '1PAG': 0.0, '1PACRF01_AE': 0.0, '1PACRF02_AG': 0.0, '1PAOther': 0.0, '1PAAll': 0.0, '1PHA': 0.0, '1PHB': 0.0, '1PHC': 0.0, '1PHD': 0.0, '1PHF': 0.0, '1PHG': 0.0, '1PHCRF01_AE': 0.0, '1PHCRF02_AG': 0.0, '1PHOther': 0.0, '1PHAll': 0.0, '2QHA': 0.0, '2QHB': 0.1, '2QHC': 0.1, '2QHD': 0.6, '2QHF': 0.0, '2QHG': 0.0, '2QHCRF01_AE': 0.4, '2QHCRF02_AG': 0.0, '2QHOther': 0.1, '2QHAll': 0.1, '2QPA': 0.0, '2QPB': 0.0, '2QPC': 0.0, '2QPD': 0.6, '2QPF': 0.0, '2QPG': 0.0, '2QPCRF01_AE': 0.0, '2QPCRF02_AG': 0.4, '2QPOther': 0.1, '2QPAll': 0.0, '2QKA': 0.0, '2QKB': 0.0, '2QKC': 0.0, '2QKD': 0.0, '2QKF': 0.0, '2QKG': 0.0, '2QKCRF01_AE': 0.0, '2QKCRF02_AG': 0.0, '2QKOther': 0.0, '2QKAll': 0.0, '2QLA': 0.0, '2QLB': 0.0, '2QLC': 0.0, '2QLD': 0.0, '2QLF': 0.2, '2QLG': 0.0, '2QLCRF01_AE': 0.0, '2QLCRF02_AG': 0.0, '2QLOther': 0.0, '2QLAll': 0.0, '2QRA': 0.0, '2QRB': 0.0, '2QRC': 0.1, '2QRD': 0.0, '2QRF': 0.0, '2QRG': 0.0, '2QRCRF01_AE': 0.0, '2QRCRF02_AG': 0.0, '2QROther': 0.0, '2QRAll': 0.0, '2QTA': 0.0, '2QTB': 0.0, '2QTC': 0.0, '2QTD': 0.0, '2QTF': 0.0, '2QTG': 0.0, '2QTCRF01_AE': 0.0, '2QTCRF02_AG': 0.0, '2QTOther': 0.0, '2QTAll': 0.0, '2QEA': 0.0, '2QEB': 0.0, '2QEC': 0.0, '2QED': 0.0, '2QEF': 0.0, '2QEG': 0.0, '2QECRF01_AE': 0.0, '2QECRF02_AG': 0.0, '2QEOther': 0.0, '2QEAll': 0.0, '2QSA': 0.0, '2QSB': 0.0, '2QSC': 0.0, '2QSD': 0.0, '2QSF': 0.0, '2QSG': 0.0, '2QSCRF01_AE': 0.0, '2QSCRF02_AG': 0.0, '2QSOther': 0.0, '2QSAll': 0.0, '2QVA': 0.0, '2QVB': 0.0, '2QVC': 0.0, '2QVD': 0.0, '2QVF': 0.0, '2QVG': 0.0, '2QVCRF01_AE': 0.0, '2QVCRF02_AG': 0.0, '2QVOther': 0.0, '2QVAll': 0.0, '2QWA': 0.0, '2QWB': 0.0, '2QWC': 0.0, '2QWD': 0.0, '2QWF': 0.0, '2QWG': 0.0, '2QWCRF01_AE': 0.0, '2QWCRF02_AG': 0.0, '2QWOther': 0.0, '2QWAll': 0.0, '2QIA': 0.0, '2QIB': 0.0, '2QIC': 0.0, '2QID': 0.0, '2QIF': 0.0, '2QIG': 0.0, '2QICRF01_AE': 0.0, '2QICRF02_AG': 0.0, '2QIOther': 0.0, '2QIAll': 0.0, '3IVA': 0.7, '3IVB': 0.4, '3IVC': 0.3, '3IVD': 0.0, '3IVF': 0.2, '3IVG': 0.0, '3IVCRF01_AE': 0.0, '3IVCRF02_AG': 0.4, '3IVOther': 0.3, '3IVAll': 0.4, '3ILA': 0.0, '3ILB': 0.0, '3ILC': 0.1, '3ILD': 0.0, '3ILF': 0.0, '3ILG': 0.3, '3ILCRF01_AE': 0.0, '3ILCRF02_AG': 0.0, '3ILOther': 0.1, '3ILAll': 0.0, '3INA': 0.0, '3INB': 0.1, '3INC': 0.1, '3IND': 0.0, '3INF': 0.2, '3ING': 0.0, '3INCRF01_AE': 0.0, '3INCRF02_AG': 0.0, '3INOther': 0.1, '3INAll': 0.1, '3IFA': 0.0, '3IFB': 0.0, '3IFC': 0.0, '3IFD': 0.0, '3IFF': 0.2, '3IFG': 0.0, '3IFCRF01_AE': 0.0, '3IFCRF02_AG': 0.0, '3IFOther': 0.0, '3IFAll': 0.0, '3IKA': 0.0, '3IKB': 0.0, '3IKC': 0.0, '3IKD': 0.0, '3IKF': 0.0, '3IKG': 0.0, '3IKCRF01_AE': 0.0, '3IKCRF02_AG': 0.0, '3IKOther': 0.0, '3IKAll': 0.0, '3ISA': 0.0, '3ISB': 0.0, '3ISC': 0.0, '3ISD': 0.0, '3ISF': 0.0, '3ISG': 0.0, '3ISCRF01_AE': 0.0, '3ISCRF02_AG': 0.0, '3ISOther': 0.0, '3ISAll': 0.0, '3ITA': 0.0, '3ITB': 0.0, '3ITC': 0.0, '3ITD': 0.0, '3ITF': 0.0, '3ITG': 0.0, '3ITCRF01_AE': 0.0, '3ITCRF02_AG': 0.0, '3ITOther': 0.0, '3ITAll': 0.0, '3IQA': 0.0, '3IQB': 0.0, '3IQC': 0.0, '3IQD': 0.0, '3IQF': 0.0, '3IQG': 0.0, '3IQCRF01_AE': 0.0, '3IQCRF02_AG': 0.0, '3IQOther': 0.0, '3IQAll': 0.0, '3IMA': 0.0, '3IMB': 0.0, '3IMC': 0.0, '3IMD': 0.0, '3IMF': 0.0, '3IMG': 0.0, '3IMCRF01_AE': 0.0, '3IMCRF02_AG': 0.0, '3IMOther': 0.0, '3IMAll': 0.0, '3IDA': 0.0, '3IDB': 0.0, '3IDC': 0.0, '3IDD': 0.0, '3IDF': 0.0, '3IDG': 0.0, '3IDCRF01_AE': 0.0, '3IDCRF02_AG': 0.0, '3IDOther': 0.0, '3IDAll': 0.0, '3IGA': 0.0, '3IGB': 0.0, '3IGC': 0.0, '3IGD': 0.0, '3IGF': 0.0, '3IGG': 0.0, '3IGCRF01_AE': 0.0, '3IGCRF02_AG': 0.0, '3IGOther': 0.0, '3IGAll': 0.0, '4TSA': 0.3, '4TSB': 0.6, '4TSC': 0.1, '4TSD': 0.4, '4TSF': 0.3, '4TSG': 0.5, '4TSCRF01_AE': 0.3, '4TSCRF02_AG': 0.0, '4TSOther': 0.4, '4TSAll': 0.5, '4TPA': 0.6, '4TPB': 0.3, '4TPC': 0.1, '4TPD': 0.0, '4TPF': 0.0, '4TPG': 0.0, '4TPCRF01_AE': 0.2, '4TPCRF02_AG': 0.4, '4TPOther': 0.2, '4TPAll': 0.2, '4TAA': 0.0, '4TAB': 0.2, '4TAC': 0.0, '4TAD': 0.0, '4TAF': 0.1, '4TAG': 0.3, '4TACRF01_AE': 0.0, '4TACRF02_AG': 0.0, '4TAOther': 0.3, '4TAAll': 0.1, '4TIA': 0.0, '4TIB': 0.0, '4TIC': 0.1, '4TID': 0.0, '4TIF': 0.1, '4TIG': 0.5, '4TICRF01_AE': 0.3, '4TICRF02_AG': 0.0, '4TIOther': 0.0, '4TIAll': 0.0, '4TNA': 0.0, '4TNB': 0.0, '4TNC': 0.0, '4TND': 0.0, '4TNF': 0.0, '4TNG': 0.3, '4TNCRF01_AE': 0.0, '4TNCRF02_AG': 0.0, '4TNOther': 0.1, '4TNAll': 0.0, '4THA': 0.0, '4THB': 0.0, '4THC': 0.0, '4THD': 0.0, '4THF': 0.1, '4THG': 0.0, '4THCRF01_AE': 0.0, '4THCRF02_AG': 0.0, '4THOther': 0.0, '4THAll': 0.0, '4TQA': 0.0, '4TQB': 0.0, '4TQC': 0.0, '4TQD': 0.0, '4TQF': 0.0, '4TQG': 0.0, '4TQCRF01_AE': 0.0, '4TQCRF02_AG': 0.0, '4TQOther': 0.0, '4TQAll': 0.0, '4TCA': 0.0, '4TCB': 0.0, '4TCC': 0.0, '4TCD': 0.0, '4TCF': 0.0, '4TCG': 0.0, '4TCCRF01_AE': 0.0, '4TCCRF02_AG': 0.0, '4TCOther': 0.0, '4TCAll': 0.0, '4TLA': 0.0, '4TLB': 0.0, '4TLC': 0.0, '4TLD': 0.0, '4TLF': 0.0, '4TLG': 0.0, '4TLCRF01_AE': 0.0, '4TLCRF02_AG': 0.0, '4TLOther': 0.0, '4TLAll': 0.0, '4TDA': 0.0, '4TDB': 0.0, '4TDC': 0.0, '4TDD': 0.0, '4TDF': 0.0, '4TDG': 0.0, '4TDCRF01_AE': 0.0, '4TDCRF02_AG': 0.0, '4TDOther': 0.0, '4TDAll': 0.0, '5LFA': 0.6, '5LFB': 0.2, '5LFC': 0.3, '5LFD': 0.9, '5LFF': 0.3, '5LFG': 0.3, '5LFCRF01_AE': 0.5, '5LFCRF02_AG': 0.0, '5LFOther': 0.2, '5LFAll': 0.2, '5LIA': 1.2, '5LIB': 0.0, '5LIC': 0.1, '5LID': 0.0, '5LIF': 0.1, '5LIG': 0.0, '5LICRF01_AE': 0.0, '5LICRF02_AG': 0.0, '5LIOther': 0.1, '5LIAll': 0.1, '5LPA': 0.0, '5LPB': 0.1, '5LPC': 0.0, '5LPD': 0.0, '5LPF': 0.1, '5LPG': 0.0, '5LPCRF01_AE': 0.0, '5LPCRF02_AG': 0.0, '5LPOther': 0.1, '5LPAll': 0.0, '5LSA': 0.0, '5LSB': 0.0, '5LSC': 0.0, '5LSD': 0.0, '5LSF': 0.1, '5LSG': 0.0, '5LSCRF01_AE': 0.0, '5LSCRF02_AG': 0.0, '5LSOther': 0.0, '5LSAll': 0.0, '5LQA': 0.0, '5LQB': 0.0, '5LQC': 0.0, '5LQD': 0.0, '5LQF': 0.0, '5LQG': 0.0, '5LQCRF01_AE': 0.0, '5LQCRF02_AG': 0.0, '5LQOther': 0.0, '5LQAll': 0.0, '5LTA': 0.0, '5LTB': 0.0, '5LTC': 0.0, '5LTD': 0.0, '5LTF': 0.0, '5LTG': 0.0, '5LTCRF01_AE': 0.0, '5LTCRF02_AG': 0.0, '5LTOther': 0.0, '5LTAll': 0.0, '5LVA': 0.0, '5LVB': 0.0, '5LVC': 0.0, '5LVD': 0.0, '5LVF': 0.0, '5LVG': 0.0, '5LVCRF01_AE': 0.0, '5LVCRF02_AG': 0.0, '5LVOther': 0.0, '5LVAll': 0.0, '5LHA': 0.0, '5LHB': 0.0, '5LHC': 0.0, '5LHD': 0.0, '5LHF': 0.0, '5LHG': 0.0, '5LHCRF01_AE': 0.0, '5LHCRF02_AG': 0.0, '5LHOther': 0.0, '5LHAll': 0.0, '5LRA': 0.0, '5LRB': 0.0, '5LRC': 0.0, '5LRD': 0.0, '5LRF': 0.0, '5LRG': 0.0, '5LRCRF01_AE': 0.0, '5LRCRF02_AG': 0.0, '5LROther': 0.0, '5LRAll': 0.0, '6WRA': 0.6, '6WRB': 0.1, '6WRC': 0.1, '6WRD': 0.4, '6WRF': 0.1, '6WRG': 0.3, '6WRCRF01_AE': 7.7, '6WRCRF02_AG': 0.0, '6WROther': 0.1, '6WRAll': 0.3, '6WLA': 0.0, '6WLB': 0.4, '6WLC': 0.2, '6WLD': 0.9, '6WLF': 0.7, '6WLG': 0.0, '6WLCRF01_AE': 0.6, '6WLCRF02_AG': 0.0, '6WLOther': 0.2, '6WLAll': 0.3, '6WGA': 0.3, '6WGB': 0.1, '6WGC': 0.1, '6WGD': 1.3, '6WGF': 0.6, '6WGG': 0.5, '6WGCRF01_AE': 0.3, '6WGCRF02_AG': 0.0, '6WGOther': 0.1, '6WGAll': 0.1, '6WCA': 0.3, '6WCB': 0.1, '6WCC': 0.0, '6WCD': 0.0, '6WCF': 0.1, '6WCG': 0.0, '6WCCRF01_AE': 0.0, '6WCCRF02_AG': 0.0, '6WCOther': 0.3, '6WCAll': 0.1, '6WVA': 0.0, '6WVB': 0.0, '6WVC': 0.0, '6WVD': 0.0, '6WVF': 0.3, '6WVG': 0.0, '6WVCRF01_AE': 0.2, '6WVCRF02_AG': 0.0, '6WVOther': 0.0, '6WVAll': 0.0, '6WTA': 0.0, '6WTB': 0.0, '6WTC': 0.0, '6WTD': 0.0, '6WTF': 0.0, '6WTG': 0.0, '6WTCRF01_AE': 0.5, '6WTCRF02_AG': 0.0, '6WTOther': 0.0, '6WTAll': 0.0, '6WSA': 0.0, '6WSB': 0.0, '6WSC': 0.0, '6WSD': 0.0, '6WSF': 0.0, '6WSG': 0.0, '6WSCRF01_AE': 0.3, '6WSCRF02_AG': 0.0, '6WSOther': 0.1, '6WSAll': 0.0, '6WMA': 0.0, '6WMB': 0.0, '6WMC': 0.0, '6WMD': 0.0, '6WMF': 0.0, '6WMG': 0.0, '6WMCRF01_AE': 0.2, '6WMCRF02_AG': 0.0, '6WMOther': 0.0, '6WMAll': 0.0, '6WFA': 0.0, '6WFB': 0.0, '6WFC': 0.0, '6WFD': 0.0, '6WFF': 0.0, '6WFG': 0.0, '6WFCRF01_AE': 0.0, '6WFCRF02_AG': 0.0, '6WFOther': 0.0, '6WFAll': 0.0, '6WEA': 0.0, '6WEB': 0.0, '6WEC': 0.0, '6WED': 0.0, '6WEF': 0.0, '6WEG': 0.0, '6WECRF01_AE': 0.0, '6WECRF02_AG': 0.0, '6WEOther': 0.0, '6WEAll': 0.0, '6WAA': 0.0, '6WAB': 0.0, '6WAC': 0.0, '6WAD': 0.0, '6WAF': 0.0, '6WAG': 0.0, '6WACRF01_AE': 0.0, '6WACRF02_AG': 0.0, '6WAOther': 0.0, '6WAAll': 0.0, '7QEA': 0.0, '7QEB': 0.2, '7QEC': 0.2, '7QED': 2.6, '7QEF': 0.1, '7QEG': 0.0, '7QECRF01_AE': 0.0, '7QECRF02_AG': 0.0, '7QEOther': 0.1, '7QEAll': 0.2, '7QHA': 0.3, '7QHB': 0.2, '7QHC': 0.0, '7QHD': 0.0, '7QHF': 0.3, '7QHG': 0.0, '7QHCRF01_AE': 0.0, '7QHCRF02_AG': 0.0, '7QHOther': 1.1, '7QHAll': 0.2, '7QLA': 0.0, '7QLB': 0.0, '7QLC': 0.0, '7QLD': 0.9, '7QLF': 0.1, '7QLG': 0.3, '7QLCRF01_AE': 0.0, '7QLCRF02_AG': 0.0, '7QLOther': 0.1, '7QLAll': 0.0, '7QDA': 0.3, '7QDB': 0.1, '7QDC': 0.0, '7QDD': 0.9, '7QDF': 0.0, '7QDG': 0.0, '7QDCRF01_AE': 0.0, '7QDCRF02_AG': 0.0, '7QDOther': 0.1, '7QDAll': 0.0, '7QRA': 0.0, '7QRB': 0.1, '7QRC': 0.2, '7QRD': 0.0, '7QRF': 0.0, '7QRG': 0.0, '7QRCRF01_AE': 0.0, '7QRCRF02_AG': 0.0, '7QROther': 0.0, '7QRAll': 0.1, '7QKA': 0.3, '7QKB': 0.1, '7QKC': 0.1, '7QKD': 0.0, '7QKF': 0.1, '7QKG': 0.0, '7QKCRF01_AE': 0.0, '7QKCRF02_AG': 0.0, '7QKOther': 0.1, '7QKAll': 0.1, '7QPA': 0.3, '7QPB': 0.0, '7QPC': 0.0, '7QPD': 0.0, '7QPF': 0.1, '7QPG': 0.0, '7QPCRF01_AE': 0.2, '7QPCRF02_AG': 0.0, '7QPOther': 0.1, '7QPAll': 0.0, '7QYA': 0.0, '7QYB': 0.0, '7QYC': 0.0, '7QYD': 0.0, '7QYF': 0.0, '7QYG': 0.0, '7QYCRF01_AE': 0.0, '7QYCRF02_AG': 0.0, '7QYOther': 0.0, '7QYAll': 0.0, '7QAA': 0.0, '7QAB': 0.0, '7QAC': 0.0, '7QAD': 0.0, '7QAF': 0.1, '7QAG': 0.0, '7QACRF01_AE': 0.0, '7QACRF02_AG': 0.0, '7QAOther': 0.0, '7QAAll': 0.0, '7QIA': 0.0, '7QIB': 0.0, '7QIC': 0.0, '7QID': 0.0, '7QIF': 0.1, '7QIG': 0.0, '7QICRF01_AE': 0.0, '7QICRF02_AG': 0.0, '7QIOther': 0.0, '7QIAll': 0.0, '7QTA': 0.0, '7QTB': 0.0, '7QTC': 0.0, '7QTD': 0.0, '7QTF': 0.0, '7QTG': 0.0, '7QTCRF01_AE': 0.0, '7QTCRF02_AG': 0.0, '7QTOther': 0.0, '7QTAll': 0.0, '7QNA': 0.0, '7QNB': 0.0, '7QNC': 0.0, '7QND': 0.0, '7QNF': 0.0, '7QNG': 0.0, '7QNCRF01_AE': 0.0, '7QNCRF02_AG': 0.0, '7QNOther': 0.0, '7QNAll': 0.0, '8RGA': 0.3, '8RGB': 0.1, '8RGC': 0.2, '8RGD': 0.0, '8RGF': 0.1, '8RGG': 1.0, '8RGCRF01_AE': 0.2, '8RGCRF02_AG': 0.0, '8RGOther': 0.0, '8RGAll': 0.1, '8RQA': 0.0, '8RQB': 0.2, '8RQC': 0.1, '8RQD': 0.4, '8RQF': 0.1, '8RQG': 0.3, '8RQCRF01_AE': 0.2, '8RQCRF02_AG': 0.0, '8RQOther': 0.1, '8RQAll': 0.2, '8RPA': 0.3, '8RPB': 0.0, '8RPC': 0.0, '8RPD': 0.0, '8RPF': 0.2, '8RPG': 0.0, '8RPCRF01_AE': 0.0, '8RPCRF02_AG': 0.0, '8RPOther': 0.0, '8RPAll': 0.0, '8RHA': 0.3, '8RHB': 0.0, '8RHC': 0.1, '8RHD': 0.0, '8RHF': 0.3, '8RHG': 0.0, '8RHCRF01_AE': 0.0, '8RHCRF02_AG': 0.0, '8RHOther': 0.0, '8RHAll': 0.0, '8RLA': 0.0, '8RLB': 0.0, '8RLC': 0.0, '8RLD': 0.4, '8RLF': 0.0, '8RLG': 0.0, '8RLCRF01_AE': 0.0, '8RLCRF02_AG': 0.0, '8RLOther': 0.0, '8RLAll': 0.0, '8RTA': 0.6, '8RTB': 0.0, '8RTC': 0.0, '8RTD': 0.0, '8RTF': 0.0, '8RTG': 0.0, '8RTCRF01_AE': 0.0, '8RTCRF02_AG': 0.0, '8RTOther': 0.0, '8RTAll': 0.0, '8RAA': 0.0, '8RAB': 0.0, '8RAC': 0.0, '8RAD': 0.0, '8RAF': 0.1, '8RAG': 0.0, '8RACRF01_AE': 0.0, '8RACRF02_AG': 0.0, '8RAOther': 0.0, '8RAAll': 0.0, '8RSA': 0.0, '8RSB': 0.0, '8RSC': 0.0, '8RSD': 0.0, '8RSF': 0.0, '8RSG': 0.0, '8RSCRF01_AE': 0.0, '8RSCRF02_AG': 0.0, '8RSOther': 0.0, '8RSAll': 0.0, '8RKA': 0.0, '8RKB': 0.0, '8RKC': 0.0, '8RKD': 0.0, '8RKF': 0.0, '8RKG': 0.0, '8RKCRF01_AE': 0.0, '8RKCRF02_AG': 0.0, '8RKOther': 0.0, '8RKAll': 0.0, '8REA': 0.0, '8REB': 0.0, '8REC': 0.0, '8RED': 0.0, '8REF': 0.0, '8REG': 0.0, '8RECRF01_AE': 0.0, '8RECRF02_AG': 0.0, '8REOther': 0.0, '8REAll': 0.0, '8R~A': 0.0, '8R~B': 0.0, '8R~C': 0.0, '8R~D': 0.0, '8R~F': 0.0, '8R~G': 0.0, '8R~CRF01_AE': 0.0, '8R~CRF02_AG': 0.0, '8R~Other': 0.0, '8R~All': 0.0, '8RWA': 0.0, '8RWB': 0.0, '8RWC': 0.0, '8RWD': 0.0, '8RWF': 0.0, '8RWG': 0.0, '8RWCRF01_AE': 0.0, '8RWCRF02_AG': 0.0, '8RWOther': 0.0, '8RWAll': 0.0, '8RDA': 0.0, '8RDB': 0.0, '8RDC': 0.0, '8RDD': 0.0, '8RDF': 0.0, '8RDG': 0.0, '8RDCRF01_AE': 0.0, '8RDCRF02_AG': 0.0, '8RDOther': 0.0, '8RDAll': 0.0, '9PTA': 0.3, '9PTB': 0.0, '9PTC': 0.0, '9PTD': 0.0, '9PTF': 0.0, '9PTG': 0.0, '9PTCRF01_AE': 0.0, '9PTCRF02_AG': 0.0, '9PTOther': 0.1, '9PTAll': 0.0, '9PLA': 0.3, '9PLB': 0.1, '9PLC': 0.1, '9PLD': 0.0, '9PLF': 0.0, '9PLG': 0.0, '9PLCRF01_AE': 0.0, '9PLCRF02_AG': 0.0, '9PLOther': 0.0, '9PLAll': 0.0, '9PAA': 0.0, '9PAB': 0.0, '9PAC': 0.0, '9PAD': 0.0, '9PAF': 0.0, '9PAG': 0.0, '9PACRF01_AE': 0.2, '9PACRF02_AG': 0.0, '9PAOther': 0.0, '9PAAll': 0.0, '9PHA': 0.3, '9PHB': 0.0, '9PHC': 0.0, '9PHD': 0.0, '9PHF': 0.0, '9PHG': 0.0, '9PHCRF01_AE': 0.0, '9PHCRF02_AG': 0.0, '9PHOther': 0.0, '9PHAll': 0.0, '9PSA': 0.0, '9PSB': 0.0, '9PSC': 0.0, '9PSD': 0.0, '9PSF': 0.1, '9PSG': 0.0, '9PSCRF01_AE': 0.0, '9PSCRF02_AG': 0.0, '9PSOther': 0.1, '9PSAll': 0.0, '9PRA': 0.0, '9PRB': 0.0, '9PRC': 0.1, '9PRD': 0.0, '9PRF': 0.0, '9PRG': 0.0, '9PRCRF01_AE': 0.0, '9PRCRF02_AG': 0.0, '9PROther': 0.1, '9PRAll': 0.0, '9PQA': 0.0, '9PQB': 0.0, '9PQC': 0.0, '9PQD': 0.0, '9PQF': 0.1, '9PQG': 0.0, '9PQCRF01_AE': 0.0, '9PQCRF02_AG': 0.0, '9PQOther': 0.0, '9PQAll': 0.0, '9PGA': 0.0, '9PGB': 0.0, '9PGC': 0.0, '9PGD': 0.0, '9PGF': 0.0, '9PGG': 0.0, '9PGCRF01_AE': 0.0, '9PGCRF02_AG': 0.0, '9PGOther': 0.0, '9PGAll': 0.0, '9PEA': 0.0, '9PEB': 0.0, '9PEC': 0.0, '9PED': 0.0, '9PEF': 0.0, '9PEG': 0.0, '9PECRF01_AE': 0.0, '9PECRF02_AG': 0.0, '9PEOther': 0.0, '9PEAll': 0.0, '9PCA': 0.0, '9PCB': 0.0, '9PCC': 0.0, '9PCD': 0.0, '9PCF': 0.0, '9PCG': 0.0, '9PCCRF01_AE': 0.0, '9PCCRF02_AG': 0.0, '9PCOther': 0.0, '9PCAll': 0.0, '10LIA': 21.0, '10LIB': 32.1, '10LIC': 6.4, '10LID': 19.4, '10LIF': 27.1, '10LIG': 22.5, '10LICRF01_AE': 21.5, '10LICRF02_AG': 16.2, '10LIOther': 28.7, '10LIAll': 27.6, '10LVA': 12.3, '10LVB': 8.2, '10LVC': 4.0, '10LVD': 17.2, '10LVF': 24.4, '10LVG': 4.7, '10LVCRF01_AE': 10.5, '10LVCRF02_AG': 20.9, '10LVOther': 23.6, '10LVAll': 10.3, '10LFA': 6.8, '10LFB': 9.5, '10LFC': 6.9, '10LFD': 7.8, '10LFF': 1.4, '10LFG': 0.5, '10LFCRF01_AE': 5.4, '10LFCRF02_AG': 3.2, '10LFOther': 0.7, '10LFAll': 7.6, '10LMA': 0.6, '10LMB': 0.0, '10LMC': 0.8, '10LMD': 0.0, '10LMF': 0.4, '10LMG': 0.8, '10LMCRF01_AE': 0.3, '10LMCRF02_AG': 0.4, '10LMOther': 0.1, '10LMAll': 0.2, '10LYA': 0.0, '10LYB': 0.4, '10LYC': 0.4, '10LYD': 0.9, '10LYF': 0.3, '10LYG': 0.0, '10LYCRF01_AE': 0.3, '10LYCRF02_AG': 0.4, '10LYOther': 0.1, '10LYAll': 0.3, '10LRA': 0.3, '10LRB': 0.5, '10LRC': 0.5, '10LRD': 0.4, '10LRF': 0.2, '10LRG': 0.0, '10LRCRF01_AE': 0.2, '10LRCRF02_AG': 0.0, '10LROther': 0.1, '10LRAll': 0.4, '10LPA': 0.6, '10LPB': 0.0, '10LPC': 0.1, '10LPD': 0.0, '10LPF': 0.4, '10LPG': 0.0, '10LPCRF01_AE': 0.2, '10LPCRF02_AG': 0.0, '10LPOther': 0.1, '10LPAll': 0.1, '10LHA': 0.3, '10LHB': 0.1, '10LHC': 0.1, '10LHD': 0.0, '10LHF': 0.1, '10LHG': 0.0, '10LHCRF01_AE': 0.2, '10LHCRF02_AG': 0.0, '10LHOther': 0.0, '10LHAll': 0.1, '10LQA': 0.0, '10LQB': 0.0, '10LQC': 0.0, '10LQD': 0.0, '10LQF': 0.1, '10LQG': 0.0, '10LQCRF01_AE': 0.2, '10LQCRF02_AG': 0.0, '10LQOther': 0.1, '10LQAll': 0.0, '10LNA': 0.0, '10LNB': 0.0, '10LNC': 0.0, '10LND': 0.0, '10LNF': 0.0, '10LNG': 0.0, '10LNCRF01_AE': 0.2, '10LNCRF02_AG': 0.0, '10LNOther': 0.0, '10LNAll': 0.0, '10LSA': 0.0, '10LSB': 0.0, '10LSC': 0.0, '10LSD': 0.0, '10LSF': 0.1, '10LSG': 0.0, '10LSCRF01_AE': 0.0, '10LSCRF02_AG': 0.0, '10LSOther': 0.0, '10LSAll': 0.0, '10LTA': 0.0, '10LTB': 0.1, '10LTC': 0.0, '10LTD': 0.0, '10LTF': 0.0, '10LTG': 0.0, '10LTCRF01_AE': 0.0, '10LTCRF02_AG': 0.0, '10LTOther': 0.0, '10LTAll': 0.0, '10LEA': 0.0, '10LEB': 0.0, '10LEC': 0.0, '10LED': 0.0, '10LEF': 0.1, '10LEG': 0.0, '10LECRF01_AE': 0.0, '10LECRF02_AG': 0.0, '10LEOther': 0.0, '10LEAll': 0.0, '10LAA': 0.0, '10LAB': 0.0, '10LAC': 0.0, '10LAD': 0.0, '10LAF': 0.0, '10LAG': 0.0, '10LACRF01_AE': 0.0, '10LACRF02_AG': 0.0, '10LAOther': 0.1, '10LAAll': 0.0, '10LGA': 0.0, '10LGB': 0.0, '10LGC': 0.0, '10LGD': 0.0, '10LGF': 0.1, '10LGG': 0.0, '10LGCRF01_AE': 0.0, '10LGCRF02_AG': 0.0, '10LGOther': 0.0, '10LGAll': 0.0, '10LKA': 0.0, '10LKB': 0.0, '10LKC': 0.0, '10LKD': 0.0, '10LKF': 0.0, '10LKG': 0.0, '10LKCRF01_AE': 0.0, '10LKCRF02_AG': 0.0, '10LKOther': 0.0, '10LKAll': 0.0, '10LCA': 0.0, '10LCB': 0.0, '10LCC': 0.0, '10LCD': 0.0, '10LCF': 0.0, '10LCG': 0.0, '10LCCRF01_AE': 0.0, '10LCCRF02_AG': 0.0, '10LCOther': 0.0, '10LCAll': 0.0, '10LDA': 0.0, '10LDB': 0.0, '10LDC': 0.0, '10LDD': 0.0, '10LDF': 0.0, '10LDG': 0.0, '10LDCRF01_AE': 0.0, '10LDCRF02_AG': 0.0, '10LDOther': 0.0, '10LDAll': 0.0, '11VIA': 0.6, '11VIB': 3.3, '11VIC': 0.7, '11VID': 1.7, '11VIF': 1.8, '11VIG': 3.7, '11VICRF01_AE': 0.4, '11VICRF02_AG': 5.1, '11VIOther': 1.4, '11VIAll': 2.7, '11VLA': 0.3, '11VLB': 0.7, '11VLC': 0.2, '11VLD': 0.4, '11VLF': 0.1, '11VLG': 0.0, '11VLCRF01_AE': 0.7, '11VLCRF02_AG': 0.4, '11VLOther': 0.2, '11VLAll': 0.6, '11VAA': 0.0, '11VAB': 0.1, '11VAC': 0.0, '11VAD': 0.0, '11VAF': 0.0, '11VAG': 0.0, '11VACRF01_AE': 0.0, '11VACRF02_AG': 0.0, '11VAOther': 0.1, '11VAAll': 0.1, '11VTA': 0.0, '11VTB': 0.0, '11VTC': 0.0, '11VTD': 0.0, '11VTF': 0.0, '11VTG': 0.0, '11VTCRF01_AE': 0.0, '11VTCRF02_AG': 0.0, '11VTOther': 0.0, '11VTAll': 0.0, '11VFA': 0.0, '11VFB': 0.0, '11VFC': 0.1, '11VFD': 0.0, '11VFF': 0.1, '11VFG': 0.0, '11VFCRF01_AE': 0.1, '11VFCRF02_AG': 0.0, '11VFOther': 0.0, '11VFAll': 0.1, '11VMA': 0.0, '11VMB': 0.0, '11VMC': 0.0, '11VMD': 0.0, '11VMF': 0.0, '11VMG': 0.0, '11VMCRF01_AE': 0.0, '11VMCRF02_AG': 0.0, '11VMOther': 0.0, '11VMAll': 0.0, '11VGA': 0.0, '11VGB': 0.0, '11VGC': 0.0, '11VGD': 0.0, '11VGF': 0.1, '11VGG': 0.0, '11VGCRF01_AE': 0.0, '11VGCRF02_AG': 0.0, '11VGOther': 0.1, '11VGAll': 0.0, '11VSA': 0.0, '11VSB': 0.0, '11VSC': 0.0, '11VSD': 0.0, '11VSF': 0.0, '11VSG': 0.0, '11VSCRF01_AE': 0.0, '11VSCRF02_AG': 0.0, '11VSOther': 0.0, '11VSAll': 0.0, '11VCA': 0.0, '11VCB': 0.0, '11VCC': 0.0, '11VCD': 0.0, '11VCF': 0.1, '11VCG': 0.0, '11VCCRF01_AE': 0.0, '11VCCRF02_AG': 0.0, '11VCOther': 0.0, '11VCAll': 0.0, '11VDA': 0.0, '11VDB': 0.0, '11VDC': 0.0, '11VDD': 0.0, '11VDF': 0.0, '11VDG': 0.0, '11VDCRF01_AE': 0.0, '11VDCRF02_AG': 0.0, '11VDOther': 0.0, '11VDAll': 0.0, '11VYA': 0.0, '11VYB': 0.0, '11VYC': 0.0, '11VYD': 0.0, '11VYF': 0.0, '11VYG': 0.0, '11VYCRF01_AE': 0.0, '11VYCRF02_AG': 0.0, '11VYOther': 0.0, '11VYAll': 0.0, '11VWA': 0.0, '11VWB': 0.0, '11VWC': 0.0, '11VWD': 0.0, '11VWF': 0.0, '11VWG': 0.0, '11VWCRF01_AE': 0.0, '11VWCRF02_AG': 0.0, '11VWOther': 0.0, '11VWAll': 0.0, '11VPA': 0.0, '11VPB': 0.0, '11VPC': 0.0, '11VPD': 0.0, '11VPF': 0.0, '11VPG': 0.0, '11VPCRF01_AE': 0.0, '11VPCRF02_AG': 0.0, '11VPOther': 0.0, '11VPAll': 0.0, '11VRA': 0.0, '11VRB': 0.0, '11VRC': 0.0, '11VRD': 0.0, '11VRF': 0.0, '11VRG': 0.0, '11VRCRF01_AE': 0.0, '11VRCRF02_AG': 0.0, '11VROther': 0.0, '11VRAll': 0.0, '12TSA': 1.2, '12TSB': 3.2, '12TSC': 44.8, '12TSD': 3.4, '12TSF': 2.7, '12TSG': 1.3, '12TSCRF01_AE': 0.0, '12TSCRF02_AG': 4.7, '12TSOther': 2.1, '12TSAll': 7.9, '12TAA': 0.9, '12TAB': 1.9, '12TAC': 1.6, '12TAD': 3.4, '12TAF': 2.6, '12TAG': 3.1, '12TACRF01_AE': 1.6, '12TACRF02_AG': 4.3, '12TAOther': 2.5, '12TAAll': 2.0, '12TPA': 1.2, '12TPB': 4.2, '12TPC': 3.1, '12TPD': 2.6, '12TPF': 1.3, '12TPG': 2.1, '12TPCRF01_AE': 1.9, '12TPCRF02_AG': 1.2, '12TPOther': 3.4, '12TPAll': 3.6, '12TKA': 0.3, '12TKB': 2.3, '12TKC': 0.1, '12TKD': 3.0, '12TKF': 3.5, '12TKG': 2.6, '12TKCRF01_AE': 0.0, '12TKCRF02_AG': 2.8, '12TKOther': 3.3, '12TKAll': 2.1, '12TIA': 0.6, '12TIB': 1.0, '12TIC': 0.1, '12TID': 1.7, '12TIF': 5.0, '12TIG': 0.5, '12TICRF01_AE': 0.1, '12TICRF02_AG': 0.0, '12TIOther': 2.8, '12TIAll': 1.2, '12TNA': 0.0, '12TNB': 2.1, '12TNC': 0.1, '12TND': 0.0, '12TNF': 3.5, '12TNG': 0.0, '12TNCRF01_AE': 0.0, '12TNCRF02_AG': 0.0, '12TNOther': 2.8, '12TNAll': 1.9, '12TEA': 0.0, '12TEB': 1.0, '12TEC': 0.1, '12TED': 0.4, '12TEF': 2.2, '12TEG': 0.3, '12TECRF01_AE': 0.0, '12TECRF02_AG': 0.4, '12TEOther': 2.2, '12TEAll': 1.0, '12TMA': 0.0, '12TMB': 0.1, '12TMC': 0.0, '12TMD': 0.4, '12TMF': 1.2, '12TMG': 0.0, '12TMCRF01_AE': 0.0, '12TMCRF02_AG': 0.4, '12TMOther': 0.9, '12TMAll': 0.2, '12TQA': 0.0, '12TQB': 0.5, '12TQC': 0.2, '12TQD': 0.0, '12TQF': 0.4, '12TQG': 0.3, '12TQCRF01_AE': 0.0, '12TQCRF02_AG': 0.0, '12TQOther': 0.5, '12TQAll': 0.4, '12TVA': 0.3, '12TVB': 0.2, '12TVC': 0.0, '12TVD': 0.0, '12TVF': 0.4, '12TVG': 0.0, '12TVCRF01_AE': 0.1, '12TVCRF02_AG': 0.8, '12TVOther': 0.4, '12TVAll': 0.2, '12TRA': 0.0, '12TRB': 0.2, '12TRC': 0.0, '12TRD': 0.0, '12TRF': 0.2, '12TRG': 0.0, '12TRCRF01_AE': 0.0, '12TRCRF02_AG': 0.0, '12TROther': 0.2, '12TRAll': 0.2, '12TDA': 0.0, '12TDB': 0.2, '12TDC': 0.0, '12TDD': 0.0, '12TDF': 0.3, '12TDG': 0.3, '12TDCRF01_AE': 0.0, '12TDCRF02_AG': 0.0, '12TDOther': 0.3, '12TDAll': 0.2, '12TLA': 0.0, '12TLB': 0.0, '12TLC': 0.0, '12TLD': 0.0, '12TLF': 0.0, '12TLG': 0.0, '12TLCRF01_AE': 0.0, '12TLCRF02_AG': 0.0, '12TLOther': 0.0, '12TLAll': 0.0, '12THA': 0.0, '12THB': 0.0, '12THC': 0.0, '12THD': 0.0, '12THF': 0.0, '12THG': 0.0, '12THCRF01_AE': 0.0, '12THCRF02_AG': 0.0, '12THOther': 0.0, '12THAll': 0.0, '12TFA': 0.0, '12TFB': 0.0, '12TFC': 0.0, '12TFD': 0.0, '12TFF': 0.0, '12TFG': 0.0, '12TFCRF01_AE': 0.0, '12TFCRF02_AG': 0.0, '12TFOther': 0.0, '12TFAll': 0.0, '12TYA': 0.0, '12TYB': 0.0, '12TYC': 0.0, '12TYD': 0.0, '12TYF': 0.0, '12TYG': 0.0, '12TYCRF01_AE': 0.0, '12TYCRF02_AG': 0.0, '12TYOther': 0.0, '12TYAll': 0.0, '12T~A': 0.0, '12T~B': 0.0, '12T~C': 0.0, '12T~D': 0.0, '12T~F': 0.0, '12T~G': 0.0, '12T~CRF01_AE': 0.0, '12T~CRF02_AG': 0.0, '12T~Other': 0.0, '12T~All': 0.0, '12TWA': 0.0, '12TWB': 0.0, '12TWC': 0.0, '12TWD': 0.0, '12TWF': 0.0, '12TWG': 0.0, '12TWCRF01_AE': 0.0, '12TWCRF02_AG': 0.0, '12TWOther': 0.0, '12TWAll': 0.0, '12TGA': 0.0, '12TGB': 0.0, '12TGC': 0.0, '12TGD': 0.0, '12TGF': 0.0, '12TGG': 0.0, '12TGCRF01_AE': 0.0, '12TGCRF02_AG': 0.0, '12TGOther': 0.0, '12TGAll': 0.0, '13IVA': 78.6, '13IVB': 28.2, '13IVC': 10.8, '13IVD': 46.8, '13IVF': 23.4, '13IVG': 95.0, '13IVCRF01_AE': 80.3, '13IVCRF02_AG': 88.1, '13IVOther': 17.2, '13IVAll': 28.9, '13IAA': 1.8, '13IAB': 0.1, '13IAC': 0.1, '13IAD': 0.0, '13IAF': 0.4, '13IAG': 4.7, '13IACRF01_AE': 0.6, '13IACRF02_AG': 10.7, '13IAOther': 0.6, '13IAAll': 0.4, '13IMA': 0.0, '13IMB': 0.5, '13IMC': 0.1, '13IMD': 0.0, '13IMF': 0.1, '13IMG': 0.0, '13IMCRF01_AE': 0.0, '13IMCRF02_AG': 0.0, '13IMOther': 0.3, '13IMAll': 0.4, '13ILA': 0.0, '13ILB': 0.2, '13ILC': 0.1, '13ILD': 0.0, '13ILF': 0.3, '13ILG': 0.0, '13ILCRF01_AE': 0.3, '13ILCRF02_AG': 0.0, '13ILOther': 0.4, '13ILAll': 0.2, '13ISA': 0.0, '13ISB': 0.0, '13ISC': 0.0, '13ISD': 0.0, '13ISF': 0.0, '13ISG': 0.0, '13ISCRF01_AE': 0.0, '13ISCRF02_AG': 0.4, '13ISOther': 0.1, '13ISAll': 0.0, '13ITA': 0.0, '13ITB': 0.0, '13ITC': 0.0, '13ITD': 0.0, '13ITF': 0.1, '13ITG': 0.0, '13ITCRF01_AE': 0.0, '13ITCRF02_AG': 0.0, '13ITOther': 0.0, '13ITAll': 0.0, '13INA': 0.0, '13INB': 0.0, '13INC': 0.0, '13IND': 0.0, '13INF': 0.1, '13ING': 0.0, '13INCRF01_AE': 0.0, '13INCRF02_AG': 0.0, '13INOther': 0.1, '13INAll': 0.0, '13IKA': 0.0, '13IKB': 0.1, '13IKC': 0.0, '13IKD': 0.0, '13IKF': 0.1, '13IKG': 0.0, '13IKCRF01_AE': 0.0, '13IKCRF02_AG': 0.0, '13IKOther': 0.0, '13IKAll': 0.0, '13IEA': 0.0, '13IEB': 0.0, '13IEC': 0.0, '13IED': 0.0, '13IEF': 0.0, '13IEG': 0.0, '13IECRF01_AE': 0.0, '13IECRF02_AG': 0.0, '13IEOther': 0.0, '13IEAll': 0.0, '13IPA': 0.0, '13IPB': 0.0, '13IPC': 0.0, '13IPD': 0.0, '13IPF': 0.0, '13IPG': 0.0, '13IPCRF01_AE': 0.0, '13IPCRF02_AG': 0.0, '13IPOther': 0.1, '13IPAll': 0.0, '13I#A': 0.0, '13I#B': 0.0, '13I#C': 0.0, '13I#D': 0.0, '13I#F': 0.0, '13I#G': 0.0, '13I#CRF01_AE': 0.0, '13I#CRF02_AG': 0.0, '13I#Other': 0.0, '13I#All': 0.0, '13IGA': 0.0, '13IGB': 0.0, '13IGC': 0.0, '13IGD': 0.0, '13IGF': 0.0, '13IGG': 0.0, '13IGCRF01_AE': 0.0, '13IGCRF02_AG': 0.0, '13IGOther': 0.0, '13IGAll': 0.0, '13IFA': 0.0, '13IFB': 0.0, '13IFC': 0.0, '13IFD': 0.0, '13IFF': 0.0, '13IFG': 0.0, '13IFCRF01_AE': 0.0, '13IFCRF02_AG': 0.0, '13IFOther': 0.0, '13IFAll': 0.0, '13IYA': 0.0, '13IYB': 0.0, '13IYC': 0.0, '13IYD': 0.0, '13IYF': 0.0, '13IYG': 0.0, '13IYCRF01_AE': 0.0, '13IYCRF02_AG': 0.0, '13IYOther': 0.0, '13IYAll': 0.0, '13I~A': 0.0, '13I~B': 0.0, '13I~C': 0.0, '13I~D': 0.0, '13I~F': 0.0, '13I~G': 0.0, '13I~CRF01_AE': 0.0, '13I~CRF02_AG': 0.0, '13I~Other': 0.0, '13I~All': 0.0, '13IQA': 0.0, '13IQB': 0.0, '13IQC': 0.0, '13IQD': 0.0, '13IQF': 0.0, '13IQG': 0.0, '13IQCRF01_AE': 0.0, '13IQCRF02_AG': 0.0, '13IQOther': 0.0, '13IQAll': 0.0, '13ICA': 0.0, '13ICB': 0.0, '13ICC': 0.0, '13ICD': 0.0, '13ICF': 0.0, '13ICG': 0.0, '13ICCRF01_AE': 0.0, '13ICCRF02_AG': 0.0, '13ICOther': 0.0, '13ICAll': 0.0, '13IRA': 0.0, '13IRB': 0.0, '13IRC': 0.0, '13IRD': 0.0, '13IRF': 0.0, '13IRG': 0.0, '13IRCRF01_AE': 0.0, '13IRCRF02_AG': 0.0, '13IROther': 0.0, '13IRAll': 0.0, '14KRA': 35.6, '14KRB': 11.0, '14KRC': 12.6, '14KRD': 8.0, '14KRF': 12.3, '14KRG': 71.1, '14KRCRF01_AE': 4.6, '14KRCRF02_AG': 59.6, '14KROther': 13.6, '14KRAll': 13.0, '14KNA': 0.3, '14KNB': 0.2, '14KNC': 0.5, '14KND': 0.0, '14KNF': 0.1, '14KNG': 0.0, '14KNCRF01_AE': 0.9, '14KNCRF02_AG': 0.4, '14KNOther': 0.1, '14KNAll': 0.3, '14KEA': 0.0, '14KEB': 0.2, '14KEC': 0.0, '14KED': 0.0, '14KEF': 0.3, '14KEG': 0.8, '14KECRF01_AE': 0.0, '14KECRF02_AG': 0.0, '14KEOther': 0.5, '14KEAll': 0.2, '14KTA': 0.0, '14KTB': 0.1, '14KTC': 0.2, '14KTD': 0.4, '14KTF': 0.5, '14KTG': 0.3, '14KTCRF01_AE': 0.1, '14KTCRF02_AG': 0.4, '14KTOther': 0.1, '14KTAll': 0.2, '14KIA': 0.0, '14KIB': 0.1, '14KIC': 0.3, '14KID': 0.0, '14KIF': 0.4, '14KIG': 0.0, '14KICRF01_AE': 0.0, '14KICRF02_AG': 0.0, '14KIOther': 1.1, '14KIAll': 0.2, '14KQA': 0.0, '14KQB': 0.0, '14KQC': 0.1, '14KQD': 0.0, '14KQF': 0.9, '14KQG': 0.0, '14KQCRF01_AE': 0.3, '14KQCRF02_AG': 0.4, '14KQOther': 0.1, '14KQAll': 0.1, '14KGA': 0.0, '14KGB': 0.0, '14KGC': 0.0, '14KGD': 0.0, '14KGF': 0.1, '14KGG': 0.0, '14KGCRF01_AE': 0.0, '14KGCRF02_AG': 0.4, '14KGOther': 0.1, '14KGAll': 0.0, '14KSA': 0.0, '14KSB': 0.0, '14KSC': 0.0, '14KSD': 0.0, '14KSF': 0.1, '14KSG': 0.0, '14KSCRF01_AE': 0.0, '14KSCRF02_AG': 0.0, '14KSOther': 0.1, '14KSAll': 0.0, '14KMA': 0.0, '14KMB': 0.2, '14KMC': 0.0, '14KMD': 0.0, '14KMF': 0.1, '14KMG': 0.0, '14KMCRF01_AE': 0.0, '14KMCRF02_AG': 0.0, '14KMOther': 0.1, '14KMAll': 0.1, '14KWA': 0.3, '14KWB': 0.0, '14KWC': 0.0, '14KWD': 0.0, '14KWF': 0.0, '14KWG': 0.0, '14KWCRF01_AE': 0.0, '14KWCRF02_AG': 0.0, '14KWOther': 0.0, '14KWAll': 0.0, '14KLA': 0.0, '14KLB': 0.0, '14KLC': 0.0, '14KLD': 0.0, '14KLF': 0.0, '14KLG': 0.0, '14KLCRF01_AE': 0.0, '14KLCRF02_AG': 0.0, '14KLOther': 0.1, '14KLAll': 0.0, '14KAA': 0.0, '14KAB': 0.0, '14KAC': 0.0, '14KAD': 0.0, '14KAF': 0.0, '14KAG': 0.0, '14KACRF01_AE': 0.0, '14KACRF02_AG': 0.0, '14KAOther': 0.0, '14KAAll': 0.0, '14KYA': 0.0, '14KYB': 0.0, '14KYC': 0.0, '14KYD': 0.0, '14KYF': 0.0, '14KYG': 0.0, '14KYCRF01_AE': 0.0, '14KYCRF02_AG': 0.0, '14KYOther': 0.0, '14KYAll': 0.0, '14KVA': 0.0, '14KVB': 0.0, '14KVC': 0.0, '14KVD': 0.0, '14KVF': 0.0, '14KVG': 0.0, '14KVCRF01_AE': 0.0, '14KVCRF02_AG': 0.0, '14KVOther': 0.0, '14KVAll': 0.0, '14KHA': 0.0, '14KHB': 0.0, '14KHC': 0.0, '14KHD': 0.0, '14KHF': 0.0, '14KHG': 0.0, '14KHCRF01_AE': 0.0, '14KHCRF02_AG': 0.0, '14KHOther': 0.0, '14KHAll': 0.0, '14KDA': 0.0, '14KDB': 0.0, '14KDC': 0.0, '14KDD': 0.0, '14KDF': 0.0, '14KDG': 0.0, '14KDCRF01_AE': 0.0, '14KDCRF02_AG': 0.0, '14KDOther': 0.0, '14KDAll': 0.0, '15IVA': 17.6, '15IVB': 25.6, '15IVC': 82.2, '15IVD': 33.1, '15IVF': 66.2, '15IVG': 20.6, '15IVCRF01_AE': 12.3, '15IVCRF02_AG': 19.8, '15IVOther': 79.2, '15IVAll': 38.2, '15ILA': 0.3, '15ILB': 0.3, '15ILC': 0.1, '15ILD': 0.8, '15ILF': 0.6, '15ILG': 0.0, '15ILCRF01_AE': 0.0, '15ILCRF02_AG': 3.5, '15ILOther': 0.1, '15ILAll': 0.3, '15IMA': 0.6, '15IMB': 0.0, '15IMC': 0.0, '15IMD': 0.0, '15IMF': 0.0, '15IMG': 0.0, '15IMCRF01_AE': 0.1, '15IMCRF02_AG': 0.0, '15IMOther': 0.0, '15IMAll': 0.0, '15ISA': 0.0, '15ISB': 0.0, '15ISC': 0.2, '15ISD': 0.0, '15ISF': 0.2, '15ISG': 0.3, '15ISCRF01_AE': 0.0, '15ISCRF02_AG': 0.0, '15ISOther': 0.1, '15ISAll': 0.1, '15INA': 0.0, '15INB': 0.0, '15INC': 0.0, '15IND': 0.0, '15INF': 0.0, '15ING': 0.0, '15INCRF01_AE': 0.0, '15INCRF02_AG': 0.4, '15INOther': 0.0, '15INAll': 0.0, '15IRA': 0.0, '15IRB': 0.0, '15IRC': 0.0, '15IRD': 0.4, '15IRF': 0.0, '15IRG': 0.0, '15IRCRF01_AE': 0.0, '15IRCRF02_AG': 0.0, '15IROther': 0.0, '15IRAll': 0.0, '15IKA': 0.0, '15IKB': 0.0, '15IKC': 0.0, '15IKD': 0.0, '15IKF': 0.1, '15IKG': 0.0, '15IKCRF01_AE': 0.0, '15IKCRF02_AG': 0.0, '15IKOther': 0.0, '15IKAll': 0.0, '15IGA': 0.0, '15IGB': 0.0, '15IGC': 0.0, '15IGD': 0.0, '15IGF': 0.1, '15IGG': 0.0, '15IGCRF01_AE': 0.0, '15IGCRF02_AG': 0.0, '15IGOther': 0.1, '15IGAll': 0.0, '15ITA': 0.0, '15ITB': 0.0, '15ITC': 0.0, '15ITD': 0.0, '15ITF': 0.0, '15ITG': 0.0, '15ITCRF01_AE': 0.0, '15ITCRF02_AG': 0.0, '15ITOther': 0.0, '15ITAll': 0.0, '15IEA': 0.0, '15IEB': 0.0, '15IEC': 0.0, '15IED': 0.0, '15IEF': 0.0, '15IEG': 0.0, '15IECRF01_AE': 0.0, '15IECRF02_AG': 0.0, '15IEOther': 0.0, '15IEAll': 0.0, '15IAA': 0.0, '15IAB': 0.0, '15IAC': 0.0, '15IAD': 0.0, '15IAF': 0.1, '15IAG': 0.0, '15IACRF01_AE': 0.0, '15IACRF02_AG': 0.0, '15IAOther': 0.0, '15IAAll': 0.0, '15IFA': 0.0, '15IFB': 0.0, '15IFC': 0.0, '15IFD': 0.0, '15IFF': 0.0, '15IFG': 0.0, '15IFCRF01_AE': 0.0, '15IFCRF02_AG': 0.0, '15IFOther': 0.0, '15IFAll': 0.0, '15IYA': 0.0, '15IYB': 0.0, '15IYC': 0.0, '15IYD': 0.0, '15IYF': 0.0, '15IYG': 0.0, '15IYCRF01_AE': 0.0, '15IYCRF02_AG': 0.0, '15IYOther': 0.0, '15IYAll': 0.0, '15ICA': 0.0, '15ICB': 0.0, '15ICC': 0.0, '15ICD': 0.0, '15ICF': 0.0, '15ICG': 0.0, '15ICCRF01_AE': 0.0, '15ICCRF02_AG': 0.0, '15ICOther': 0.0, '15ICAll': 0.0, '16GEA': 26.9, '16GEB': 4.8, '16GEC': 10.9, '16GED': 7.1, '16GEF': 17.7, '16GEG': 4.4, '16GECRF01_AE': 17.2, '16GECRF02_AG': 30.8, '16GEOther': 11.8, '16GEAll': 7.7, '16GAA': 3.3, '16GAB': 3.0, '16GAC': 2.8, '16GAD': 0.8, '16GAF': 1.4, '16GAG': 1.3, '16GACRF01_AE': 2.5, '16GACRF02_AG': 0.4, '16GAOther': 2.1, '16GAAll': 2.7, '16GRA': 0.0, '16GRB': 0.2, '16GRC': 0.2, '16GRD': 0.0, '16GRF': 0.3, '16GRG': 0.3, '16GRCRF01_AE': 0.1, '16GRCRF02_AG': 0.4, '16GROther': 0.2, '16GRAll': 0.2, '16GQA': 0.0, '16GQB': 0.1, '16GQC': 0.0, '16GQD': 0.4, '16GQF': 0.6, '16GQG': 0.0, '16GQCRF01_AE': 0.1, '16GQCRF02_AG': 0.4, '16GQOther': 0.2, '16GQAll': 0.1, '16GWA': 0.0, '16GWB': 0.0, '16GWC': 0.1, '16GWD': 0.0, '16GWF': 0.1, '16GWG': 0.0, '16GWCRF01_AE': 0.0, '16GWCRF02_AG': 0.0, '16GWOther': 0.1, '16GWAll': 0.0, '16GDA': 0.0, '16GDB': 0.0, '16GDC': 0.0, '16GDD': 0.0, '16GDF': 0.0, '16GDG': 0.0, '16GDCRF01_AE': 0.0, '16GDCRF02_AG': 0.0, '16GDOther': 0.0, '16GDAll': 0.0, '16GSA': 0.0, '16GSB': 0.0, '16GSC': 0.0, '16GSD': 0.0, '16GSF': 0.0, '16GSG': 0.0, '16GSCRF01_AE': 0.0, '16GSCRF02_AG': 0.0, '16GSOther': 0.0, '16GSAll': 0.0, '16GKA': 0.0, '16GKB': 0.0, '16GKC': 0.0, '16GKD': 0.0, '16GKF': 0.0, '16GKG': 0.0, '16GKCRF01_AE': 0.0, '16GKCRF02_AG': 0.0, '16GKOther': 0.0, '16GKAll': 0.0, '16GYA': 0.0, '16GYB': 0.0, '16GYC': 0.0, '16GYD': 0.0, '16GYF': 0.0, '16GYG': 0.0, '16GYCRF01_AE': 0.0, '16GYCRF02_AG': 0.0, '16GYOther': 0.0, '16GYAll': 0.0, '16GVA': 0.0, '16GVB': 0.0, '16GVC': 0.0, '16GVD': 0.0, '16GVF': 0.0, '16GVG': 0.0, '16GVCRF01_AE': 0.0, '16GVCRF02_AG': 0.0, '16GVOther': 0.0, '16GVAll': 0.0, '16GMA': 0.0, '16GMB': 0.0, '16GMC': 0.0, '16GMD': 0.0, '16GMF': 0.0, '16GMG': 0.0, '16GMCRF01_AE': 0.0, '16GMCRF02_AG': 0.0, '16GMOther': 0.0, '16GMAll': 0.0, '16GCA': 0.0, '16GCB': 0.0, '16GCC': 0.0, '16GCD': 0.0, '16GCF': 0.0, '16GCG': 0.0, '16GCCRF01_AE': 0.0, '16GCCRF02_AG': 0.0, '16GCOther': 0.0, '16GCAll': 0.0, '16GLA': 0.0, '16GLB': 0.0, '16GLC': 0.0, '16GLD': 0.0, '16GLF': 0.0, '16GLG': 0.0, '16GLCRF01_AE': 0.0, '16GLCRF02_AG': 0.0, '16GLOther': 0.0, '16GLAll': 0.0, '17GEA': 0.6, '17GEB': 1.5, '17GEC': 0.0, '17GED': 0.0, '17GEF': 3.4, '17GEG': 3.9, '17GECRF01_AE': 1.5, '17GECRF02_AG': 6.6, '17GEOther': 3.3, '17GEAll': 1.6, '17GDA': 0.6, '17GDB': 1.1, '17GDC': 1.1, '17GDD': 0.8, '17GDF': 1.4, '17GDG': 0.0, '17GDCRF01_AE': 0.3, '17GDCRF02_AG': 0.0, '17GDOther': 2.2, '17GDAll': 1.1, '17GRA': 0.3, '17GRB': 0.1, '17GRC': 0.0, '17GRD': 0.0, '17GRF': 0.1, '17GRG': 0.0, '17GRCRF01_AE': 0.1, '17GRCRF02_AG': 0.0, '17GROther': 0.1, '17GRAll': 0.1, '17GAA': 0.0, '17GAB': 0.0, '17GAC': 0.0, '17GAD': 0.0, '17GAF': 0.1, '17GAG': 0.0, '17GACRF01_AE': 0.1, '17GACRF02_AG': 0.0, '17GAOther': 0.0, '17GAAll': 0.0, '17GSA': 0.0, '17GSB': 0.0, '17GSC': 0.1, '17GSD': 0.0, '17GSF': 0.0, '17GSG': 0.0, '17GSCRF01_AE': 0.0, '17GSCRF02_AG': 0.0, '17GSOther': 0.0, '17GSAll': 0.0, '17GKA': 0.0, '17GKB': 0.0, '17GKC': 0.0, '17GKD': 0.0, '17GKF': 0.0, '17GKG': 0.0, '17GKCRF01_AE': 0.0, '17GKCRF02_AG': 0.0, '17GKOther': 0.0, '17GKAll': 0.0, '17GNA': 0.0, '17GNB': 0.0, '17GNC': 0.0, '17GND': 0.0, '17GNF': 0.0, '17GNG': 0.0, '17GNCRF01_AE': 0.0, '17GNCRF02_AG': 0.0, '17GNOther': 0.0, '17GNAll': 0.0, '17GVA': 0.0, '17GVB': 0.0, '17GVC': 0.0, '17GVD': 0.0, '17GVF': 0.0, '17GVG': 0.0, '17GVCRF01_AE': 0.1, '17GVCRF02_AG': 0.0, '17GVOther': 0.0, '17GVAll': 0.0, '17GCA': 0.0, '17GCB': 0.0, '17GCC': 0.0, '17GCD': 0.0, '17GCF': 0.0, '17GCG': 0.0, '17GCCRF01_AE': 0.0, '17GCCRF02_AG': 0.0, '17GCOther': 0.0, '17GCAll': 0.0, '17GWA': 0.0, '17GWB': 0.0, '17GWC': 0.0, '17GWD': 0.0, '17GWF': 0.0, '17GWG': 0.0, '17GWCRF01_AE': 0.0, '17GWCRF02_AG': 0.0, '17GWOther': 0.0, '17GWAll': 0.0, '18QEA': 0.3, '18QEB': 0.5, '18QEC': 0.8, '18QED': 0.4, '18QEF': 1.8, '18QEG': 0.0, '18QECRF01_AE': 0.4, '18QECRF02_AG': 0.0, '18QEOther': 1.6, '18QEAll': 0.7, '18QHA': 0.3, '18QHB': 2.0, '18QHC': 0.1, '18QHD': 0.4, '18QHF': 0.6, '18QHG': 0.0, '18QHCRF01_AE': 0.3, '18QHCRF02_AG': 0.0, '18QHOther': 0.7, '18QHAll': 1.4, '18QLA': 0.0, '18QLB': 0.4, '18QLC': 0.2, '18QLD': 1.3, '18QLF': 1.1, '18QLG': 0.0, '18QLCRF01_AE': 0.0, '18QLCRF02_AG': 0.0, '18QLOther': 0.8, '18QLAll': 0.4, '18QKA': 0.0, '18QKB': 0.2, '18QKC': 0.0, '18QKD': 0.4, '18QKF': 0.3, '18QKG': 0.3, '18QKCRF01_AE': 0.0, '18QKCRF02_AG': 0.8, '18QKOther': 0.1, '18QKAll': 0.2, '18QRA': 0.3, '18QRB': 0.2, '18QRC': 0.1, '18QRD': 0.4, '18QRF': 0.1, '18QRG': 0.0, '18QRCRF01_AE': 0.3, '18QRCRF02_AG': 0.0, '18QROther': 0.1, '18QRAll': 0.1, '18QIA': 0.0, '18QIB': 0.2, '18QIC': 0.0, '18QID': 0.0, '18QIF': 0.0, '18QIG': 0.0, '18QICRF01_AE': 0.0, '18QICRF02_AG': 0.4, '18QIOther': 0.1, '18QIAll': 0.1, '18QMA': 0.0, '18QMB': 0.0, '18QMC': 0.0, '18QMD': 0.0, '18QMF': 0.3, '18QMG': 0.0, '18QMCRF01_AE': 0.0, '18QMCRF02_AG': 0.0, '18QMOther': 0.1, '18QMAll': 0.0, '18QSA': 0.3, '18QSB': 0.0, '18QSC': 0.0, '18QSD': 0.0, '18QSF': 0.0, '18QSG': 0.0, '18QSCRF01_AE': 0.0, '18QSCRF02_AG': 0.0, '18QSOther': 0.0, '18QSAll': 0.0, '18QPA': 0.0, '18QPB': 0.0, '18QPC': 0.1, '18QPD': 0.0, '18QPF': 0.0, '18QPG': 0.0, '18QPCRF01_AE': 0.0, '18QPCRF02_AG': 0.0, '18QPOther': 0.0, '18QPAll': 0.0, '18QTA': 0.0, '18QTB': 0.0, '18QTC': 0.0, '18QTD': 0.0, '18QTF': 0.0, '18QTG': 0.0, '18QTCRF01_AE': 0.0, '18QTCRF02_AG': 0.0, '18QTOther': 0.1, '18QTAll': 0.0, '18QYA': 0.0, '18QYB': 0.0, '18QYC': 0.0, '18QYD': 0.0, '18QYF': 0.1, '18QYG': 0.0, '18QYCRF01_AE': 0.0, '18QYCRF02_AG': 0.0, '18QYOther': 0.0, '18QYAll': 0.0, '18QVA': 0.0, '18QVB': 0.0, '18QVC': 0.0, '18QVD': 0.0, '18QVF': 0.1, '18QVG': 0.0, '18QVCRF01_AE': 0.0, '18QVCRF02_AG': 0.0, '18QVOther': 0.0, '18QVAll': 0.0, '18QAA': 0.0, '18QAB': 0.0, '18QAC': 0.0, '18QAD': 0.0, '18QAF': 0.0, '18QAG': 0.0, '18QACRF01_AE': 0.0, '18QACRF02_AG': 0.0, '18QAOther': 0.0, '18QAAll': 0.0, '18QGA': 0.0, '18QGB': 0.0, '18QGC': 0.0, '18QGD': 0.0, '18QGF': 0.0, '18QGG': 0.0, '18QGCRF01_AE': 0.0, '18QGCRF02_AG': 0.0, '18QGOther': 0.0, '18QGAll': 0.0, '18QNA': 0.0, '18QNB': 0.0, '18QNC': 0.0, '18QND': 0.0, '18QNF': 0.0, '18QNG': 0.0, '18QNCRF01_AE': 0.0, '18QNCRF02_AG': 0.0, '18QNOther': 0.0, '18QNAll': 0.0, '18QWA': 0.0, '18QWB': 0.0, '18QWC': 0.0, '18QWD': 0.0, '18QWF': 0.0, '18QWG': 0.0, '18QWCRF01_AE': 0.0, '18QWCRF02_AG': 0.0, '18QWOther': 0.0, '18QWAll': 0.0, '18QDA': 0.0, '18QDB': 0.0, '18QDC': 0.0, '18QDD': 0.0, '18QDF': 0.0, '18QDG': 0.0, '18QDCRF01_AE': 0.0, '18QDCRF02_AG': 0.0, '18QDOther': 0.0, '18QDAll': 0.0, '19LIA': 2.4, '19LIB': 11.8, '19LIC': 47.2, '19LID': 12.6, '19LIF': 10.1, '19LIG': 4.2, '19LICRF01_AE': 1.0, '19LICRF02_AG': 3.8, '19LIOther': 5.8, '19LIAll': 14.7, '19LPA': 0.9, '19LPB': 0.9, '19LPC': 0.3, '19LPD': 0.0, '19LPF': 0.9, '19LPG': 2.6, '19LPCRF01_AE': 0.3, '19LPCRF02_AG': 15.8, '19LPOther': 0.8, '19LPAll': 1.0, '19LTA': 0.0, '19LTB': 1.3, '19LTC': 11.8, '19LTD': 1.3, '19LTF': 0.3, '19LTG': 0.5, '19LTCRF01_AE': 0.0, '19LTCRF02_AG': 0.8, '19LTOther': 0.2, '19LTAll': 2.3, '19LVA': 0.6, '19LVB': 2.4, '19LVC': 5.9, '19LVD': 2.5, '19LVF': 2.1, '19LVG': 0.0, '19LVCRF01_AE': 0.7, '19LVCRF02_AG': 0.8, '19LVOther': 1.4, '19LVAll': 2.6, '19LQA': 0.0, '19LQB': 1.9, '19LQC': 1.0, '19LQD': 2.1, '19LQF': 1.2, '19LQG': 0.0, '19LQCRF01_AE': 0.0, '19LQCRF02_AG': 0.0, '19LQOther': 2.4, '19LQAll': 1.7, '19LMA': 0.0, '19LMB': 0.2, '19LMC': 0.5, '19LMD': 0.0, '19LMF': 0.0, '19LMG': 0.3, '19LMCRF01_AE': 0.4, '19LMCRF02_AG': 0.0, '19LMOther': 0.2, '19LMAll': 0.2, '19LSA': 0.0, '19LSB': 0.1, '19LSC': 0.0, '19LSD': 0.4, '19LSF': 0.4, '19LSG': 1.0, '19LSCRF01_AE': 0.0, '19LSCRF02_AG': 0.0, '19LSOther': 0.2, '19LSAll': 0.1, '19LKA': 0.0, '19LKB': 0.0, '19LKC': 0.2, '19LKD': 0.0, '19LKF': 0.0, '19LKG': 0.0, '19LKCRF01_AE': 0.0, '19LKCRF02_AG': 0.4, '19LKOther': 0.0, '19LKAll': 0.0, '19LRA': 0.0, '19LRB': 0.1, '19LRC': 0.0, '19LRD': 0.0, '19LRF': 0.0, '19LRG': 0.5, '19LRCRF01_AE': 0.0, '19LRCRF02_AG': 0.0, '19LROther': 0.1, '19LRAll': 0.1, '19LEA': 0.0, '19LEB': 0.1, '19LEC': 0.2, '19LED': 0.0, '19LEF': 0.0, '19LEG': 0.0, '19LECRF01_AE': 0.1, '19LECRF02_AG': 0.0, '19LEOther': 0.0, '19LEAll': 0.1, '19LFA': 0.3, '19LFB': 0.1, '19LFC': 0.0, '19LFD': 0.0, '19LFF': 0.1, '19LFG': 0.0, '19LFCRF01_AE': 0.0, '19LFCRF02_AG': 0.0, '19LFOther': 0.0, '19LFAll': 0.1, '19LHA': 0.0, '19LHB': 0.0, '19LHC': 0.0, '19LHD': 0.0, '19LHF': 0.0, '19LHG': 0.0, '19LHCRF01_AE': 0.0, '19LHCRF02_AG': 0.0, '19LHOther': 0.0, '19LHAll': 0.0, '19LNA': 0.0, '19LNB': 0.0, '19LNC': 0.0, '19LND': 0.0, '19LNF': 0.0, '19LNG': 0.0, '19LNCRF01_AE': 0.0, '19LNCRF02_AG': 0.0, '19LNOther': 0.0, '19LNAll': 0.0, '19LAA': 0.0, '19LAB': 0.0, '19LAC': 0.0, '19LAD': 0.0, '19LAF': 0.0, '19LAG': 0.0, '19LACRF01_AE': 0.0, '19LACRF02_AG': 0.0, '19LAOther': 0.0, '19LAAll': 0.0, '19LWA': 0.0, '19LWB': 0.0, '19LWC': 0.0, '19LWD': 0.0, '19LWF': 0.0, '19LWG': 0.0, '19LWCRF01_AE': 0.0, '19LWCRF02_AG': 0.0, '19LWOther': 0.0, '19LWAll': 0.0, '19LDA': 0.0, '19LDB': 0.0, '19LDC': 0.0, '19LDD': 0.0, '19LDF': 0.0, '19LDG': 0.0, '19LDCRF01_AE': 0.0, '19LDCRF02_AG': 0.0, '19LDOther': 0.0, '19LDAll': 0.0, '19LGA': 0.0, '19LGB': 0.0, '19LGC': 0.0, '19LGD': 0.0, '19LGF': 0.0, '19LGG': 0.0, '19LGCRF01_AE': 0.0, '19LGCRF02_AG': 0.0, '19LGOther': 0.0, '19LGAll': 0.0, '20KIA': 11.6, '20KIB': 6.4, '20KIC': 1.3, '20KID': 6.3, '20KIF': 2.6, '20KIG': 91.4, '20KICRF01_AE': 10.8, '20KICRF02_AG': 91.2, '20KIOther': 4.4, '20KIAll': 7.8, '20KRA': 17.6, '20KRB': 13.6, '20KRC': 28.2, '20KRD': 18.0, '20KRF': 45.0, '20KRG': 0.3, '20KRCRF01_AE': 19.7, '20KRCRF02_AG': 4.6, '20KROther': 33.5, '20KRAll': 18.7, '20KTA': 4.2, '20KTB': 4.4, '20KTC': 5.0, '20KTD': 10.9, '20KTF': 10.0, '20KTG': 3.4, '20KTCRF01_AE': 3.3, '20KTCRF02_AG': 0.0, '20KTOther': 10.3, '20KTAll': 5.3, '20KMA': 0.3, '20KMB': 2.6, '20KMC': 1.2, '20KMD': 3.8, '20KMF': 3.4, '20KMG': 0.0, '20KMCRF01_AE': 0.1, '20KMCRF02_AG': 0.4, '20KMOther': 5.2, '20KMAll': 2.5, '20KVA': 0.6, '20KVB': 0.9, '20KVC': 0.1, '20KVD': 1.3, '20KVF': 0.8, '20KVG': 4.2, '20KVCRF01_AE': 1.0, '20KVCRF02_AG': 3.1, '20KVOther': 0.8, '20KVAll': 0.8, '20KLA': 0.0, '20KLB': 0.2, '20KLC': 0.0, '20KLD': 0.0, '20KLF': 0.1, '20KLG': 0.3, '20KLCRF01_AE': 0.0, '20KLCRF02_AG': 0.0, '20KLOther': 0.1, '20KLAll': 0.2, '20KQA': 0.0, '20KQB': 0.0, '20KQC': 0.1, '20KQD': 0.0, '20KQF': 0.1, '20KQG': 0.0, '20KQCRF01_AE': 0.0, '20KQCRF02_AG': 0.0, '20KQOther': 0.1, '20KQAll': 0.1, '20KNA': 0.0, '20KNB': 0.0, '20KNC': 0.0, '20KND': 0.0, '20KNF': 0.0, '20KNG': 0.3, '20KNCRF01_AE': 0.1, '20KNCRF02_AG': 0.0, '20KNOther': 0.0, '20KNAll': 0.0, '20KAA': 0.0, '20KAB': 0.0, '20KAC': 0.0, '20KAD': 0.0, '20KAF': 0.1, '20KAG': 0.0, '20KACRF01_AE': 0.0, '20KACRF02_AG': 0.0, '20KAOther': 0.1, '20KAAll': 0.0, '20KPA': 0.0, '20KPB': 0.0, '20KPC': 0.0, '20KPD': 0.0, '20KPF': 0.1, '20KPG': 0.0, '20KPCRF01_AE': 0.0, '20KPCRF02_AG': 0.0, '20KPOther': 0.1, '20KPAll': 0.0, '20KEA': 0.0, '20KEB': 0.0, '20KEC': 0.0, '20KED': 0.0, '20KEF': 0.0, '20KEG': 0.0, '20KECRF01_AE': 0.0, '20KECRF02_AG': 0.0, '20KEOther': 0.0, '20KEAll': 0.0, '20KDA': 0.0, '20KDB': 0.0, '20KDC': 0.0, '20KDD': 0.0, '20KDF': 0.0, '20KDG': 0.0, '20KDCRF01_AE': 0.0, '20KDCRF02_AG': 0.0, '20KDOther': 0.0, '20KDAll': 0.0, '20KGA': 0.0, '20KGB': 0.0, '20KGC': 0.0, '20KGD': 0.0, '20KGF': 0.0, '20KGG': 0.0, '20KGCRF01_AE': 0.0, '20KGCRF02_AG': 0.0, '20KGOther': 0.0, '20KGAll': 0.0, '21EKA': 0.6, '21EKB': 0.2, '21EKC': 0.4, '21EKD': 0.8, '21EKF': 0.3, '21EKG': 0.8, '21EKCRF01_AE': 0.6, '21EKCRF02_AG': 0.4, '21EKOther': 0.2, '21EKAll': 0.3, '21EDA': 0.0, '21EDB': 0.1, '21EDC': 0.0, '21EDD': 0.0, '21EDF': 0.1, '21EDG': 0.0, '21EDCRF01_AE': 0.1, '21EDCRF02_AG': 1.2, '21EDOther': 0.3, '21EDAll': 0.1, '21EQA': 0.3, '21EQB': 0.1, '21EQC': 0.0, '21EQD': 0.0, '21EQF': 0.3, '21EQG': 0.3, '21EQCRF01_AE': 0.0, '21EQCRF02_AG': 0.0, '21EQOther': 0.1, '21EQAll': 0.1, '21EGA': 0.0, '21EGB': 0.0, '21EGC': 0.0, '21EGD': 0.0, '21EGF': 0.1, '21EGG': 0.0, '21EGCRF01_AE': 0.0, '21EGCRF02_AG': 0.0, '21EGOther': 0.1, '21EGAll': 0.0, '21EVA': 0.0, '21EVB': 0.1, '21EVC': 0.0, '21EVD': 0.0, '21EVF': 0.0, '21EVG': 0.0, '21EVCRF01_AE': 0.0, '21EVCRF02_AG': 0.0, '21EVOther': 0.1, '21EVAll': 0.0, '21ETA': 0.0, '21ETB': 0.0, '21ETC': 0.0, '21ETD': 0.0, '21ETF': 0.0, '21ETG': 0.0, '21ETCRF01_AE': 0.0, '21ETCRF02_AG': 0.0, '21ETOther': 0.0, '21ETAll': 0.0, '21EAA': 0.0, '21EAB': 0.0, '21EAC': 0.0, '21EAD': 0.0, '21EAF': 0.0, '21EAG': 0.0, '21EACRF01_AE': 0.0, '21EACRF02_AG': 0.0, '21EAOther': 0.1, '21EAAll': 0.0, '21ERA': 0.0, '21ERB': 0.0, '21ERC': 0.0, '21ERD': 0.0, '21ERF': 0.0, '21ERG': 0.0, '21ERCRF01_AE': 0.0, '21ERCRF02_AG': 0.0, '21EROther': 0.1, '21ERAll': 0.0, '21EIA': 0.0, '21EIB': 0.0, '21EIC': 0.0, '21EID': 0.0, '21EIF': 0.0, '21EIG': 0.0, '21EICRF01_AE': 0.0, '21EICRF02_AG': 0.0, '21EIOther': 0.0, '21EIAll': 0.0, '22AVA': 2.7, '22AVB': 0.7, '22AVC': 0.2, '22AVD': 0.8, '22AVF': 1.1, '22AVG': 2.3, '22AVCRF01_AE': 2.2, '22AVCRF02_AG': 5.8, '22AVOther': 0.5, '22AVAll': 0.8, '22ATA': 0.0, '22ATB': 0.1, '22ATC': 0.0, '22ATD': 0.0, '22ATF': 0.0, '22ATG': 0.0, '22ATCRF01_AE': 0.0, '22ATCRF02_AG': 0.4, '22ATOther': 0.1, '22ATAll': 0.1, '22ARA': 0.3, '22ARB': 0.0, '22ARC': 0.0, '22ARD': 0.0, '22ARF': 0.0, '22ARG': 0.0, '22ARCRF01_AE': 0.0, '22ARCRF02_AG': 0.0, '22AROther': 0.0, '22ARAll': 0.0, '22ASA': 0.0, '22ASB': 0.0, '22ASC': 0.0, '22ASD': 0.0, '22ASF': 0.0, '22ASG': 0.0, '22ASCRF01_AE': 0.0, '22ASCRF02_AG': 0.0, '22ASOther': 0.0, '22ASAll': 0.0, '22APA': 0.0, '22APB': 0.0, '22APC': 0.0, '22APD': 0.0, '22APF': 0.0, '22APG': 0.0, '22APCRF01_AE': 0.0, '22APCRF02_AG': 0.0, '22APOther': 0.0, '22APAll': 0.0, '22ADA': 0.0, '22ADB': 0.0, '22ADC': 0.0, '22ADD': 0.0, '22ADF': 0.0, '22ADG': 0.0, '22ADCRF01_AE': 0.0, '22ADCRF02_AG': 0.0, '22ADOther': 0.1, '22ADAll': 0.0, '22ACA': 0.0, '22ACB': 0.0, '22ACC': 0.0, '22ACD': 0.0, '22ACF': 0.0, '22ACG': 0.0, '22ACCRF01_AE': 0.0, '22ACCRF02_AG': 0.0, '22ACOther': 0.0, '22ACAll': 0.0, '22ALA': 0.0, '22ALB': 0.0, '22ALC': 0.0, '22ALD': 0.0, '22ALF': 0.0, '22ALG': 0.0, '22ALCRF01_AE': 0.0, '22ALCRF02_AG': 0.0, '22ALOther': 0.0, '22ALAll': 0.0, '22AGA': 0.0, '22AGB': 0.0, '22AGC': 0.0, '22AGD': 0.0, '22AGF': 0.0, '22AGG': 0.0, '22AGCRF01_AE': 0.0, '22AGCRF02_AG': 0.0, '22AGOther': 0.0, '22AGAll': 0.0, '23LIA': 0.6, '23LIB': 1.4, '23LIC': 1.9, '23LID': 1.7, '23LIF': 1.0, '23LIG': 1.6, '23LICRF01_AE': 0.3, '23LICRF02_AG': 0.4, '23LIOther': 0.6, '23LIAll': 1.3, '23LFA': 0.0, '23LFB': 0.0, '23LFC': 0.5, '23LFD': 0.0, '23LFF': 0.0, '23LFG': 0.0, '23LFCRF01_AE': 0.0, '23LFCRF02_AG': 0.0, '23LFOther': 0.1, '23LFAll': 0.1, '23LVA': 0.0, '23LVB': 0.1, '23LVC': 0.0, '23LVD': 0.0, '23LVF': 0.3, '23LVG': 0.0, '23LVCRF01_AE': 0.0, '23LVCRF02_AG': 0.0, '23LVOther': 0.1, '23LVAll': 0.1, '23LKA': 0.0, '23LKB': 0.0, '23LKC': 0.0, '23LKD': 0.4, '23LKF': 0.0, '23LKG': 0.0, '23LKCRF01_AE': 0.0, '23LKCRF02_AG': 0.0, '23LKOther': 0.0, '23LKAll': 0.0, '23LSA': 0.0, '23LSB': 0.0, '23LSC': 0.0, '23LSD': 0.0, '23LSF': 0.1, '23LSG': 0.0, '23LSCRF01_AE': 0.0, '23LSCRF02_AG': 0.0, '23LSOther': 0.1, '23LSAll': 0.0, '23LPA': 0.0, '23LPB': 0.0, '23LPC': 0.0, '23LPD': 0.0, '23LPF': 0.0, '23LPG': 0.0, '23LPCRF01_AE': 0.0, '23LPCRF02_AG': 0.0, '23LPOther': 0.0, '23LPAll': 0.0, '23LQA': 0.0, '23LQB': 0.0, '23LQC': 0.0, '23LQD': 0.0, '23LQF': 0.0, '23LQG': 0.0, '23LQCRF01_AE': 0.0, '23LQCRF02_AG': 0.0, '23LQOther': 0.0, '23LQAll': 0.0, '23LMA': 0.0, '23LMB': 0.0, '23LMC': 0.0, '23LMD': 0.0, '23LMF': 0.0, '23LMG': 0.0, '23LMCRF01_AE': 0.0, '23LMCRF02_AG': 0.0, '23LMOther': 0.1, '23LMAll': 0.0, '23LYA': 0.0, '23LYB': 0.0, '23LYC': 0.0, '23LYD': 0.0, '23LYF': 0.0, '23LYG': 0.0, '23LYCRF01_AE': 0.0, '23LYCRF02_AG': 0.0, '23LYOther': 0.0, '23LYAll': 0.0, '23LHA': 0.0, '23LHB': 0.0, '23LHC': 0.0, '23LHD': 0.0, '23LHF': 0.0, '23LHG': 0.0, '23LHCRF01_AE': 0.0, '23LHCRF02_AG': 0.0, '23LHOther': 0.0, '23LHAll': 0.0, '23LRA': 0.0, '23LRB': 0.0, '23LRC': 0.0, '23LRD': 0.0, '23LRF': 0.0, '23LRG': 0.0, '23LRCRF01_AE': 0.0, '23LRCRF02_AG': 0.0, '23LROther': 0.0, '23LRAll': 0.0, '24LIA': 1.5, '24LIB': 5.6, '24LIC': 2.2, '24LID': 5.9, '24LIF': 8.5, '24LIG': 2.9, '24LICRF01_AE': 0.4, '24LICRF02_AG': 0.4, '24LIOther': 5.3, '24LIAll': 5.1, '24LFA': 0.0, '24LFB': 0.7, '24LFC': 0.1, '24LFD': 0.0, '24LFF': 0.1, '24LFG': 0.0, '24LFCRF01_AE': 0.0, '24LFCRF02_AG': 0.0, '24LFOther': 0.4, '24LFAll': 0.5, '24LSA': 0.0, '24LSB': 0.1, '24LSC': 0.0, '24LSD': 0.0, '24LSF': 0.0, '24LSG': 0.0, '24LSCRF01_AE': 0.1, '24LSCRF02_AG': 0.0, '24LSOther': 0.0, '24LSAll': 0.0, '24LMA': 0.0, '24LMB': 0.1, '24LMC': 0.0, '24LMD': 0.0, '24LMF': 0.0, '24LMG': 0.0, '24LMCRF01_AE': 0.0, '24LMCRF02_AG': 0.0, '24LMOther': 0.1, '24LMAll': 0.1, '24LVA': 0.3, '24LVB': 0.0, '24LVC': 0.0, '24LVD': 0.0, '24LVF': 0.0, '24LVG': 0.0, '24LVCRF01_AE': 0.0, '24LVCRF02_AG': 0.0, '24LVOther': 0.0, '24LVAll': 0.0, '24LKA': 0.0, '24LKB': 0.0, '24LKC': 0.0, '24LKD': 0.0, '24LKF': 0.1, '24LKG': 0.0, '24LKCRF01_AE': 0.0, '24LKCRF02_AG': 0.0, '24LKOther': 0.0, '24LKAll': 0.0, '24L~A': 0.0, '24L~B': 0.0, '24L~C': 0.0, '24L~D': 0.0, '24L~F': 0.1, '24L~G': 0.0, '24L~CRF01_AE': 0.0, '24L~CRF02_AG': 0.0, '24L~Other': 0.0, '24L~All': 0.0, '24LRA': 0.0, '24LRB': 0.0, '24LRC': 0.0, '24LRD': 0.0, '24LRF': 0.0, '24LRG': 0.0, '24LRCRF01_AE': 0.0, '24LRCRF02_AG': 0.0, '24LROther': 0.1, '24LRAll': 0.0, '24LWA': 0.0, '24LWB': 0.0, '24LWC': 0.0, '24LWD': 0.0, '24LWF': 0.0, '24LWG': 0.0, '24LWCRF01_AE': 0.0, '24LWCRF02_AG': 0.0, '24LWOther': 0.0, '24LWAll': 0.0, '24LPA': 0.0, '24LPB': 0.0, '24LPC': 0.0, '24LPD': 0.0, '24LPF': 0.0, '24LPG': 0.0, '24LPCRF01_AE': 0.0, '24LPCRF02_AG': 0.0, '24LPOther': 0.0, '24LPAll': 0.0, '24LGA': 0.0, '24LGB': 0.0, '24LGC': 0.0, '24LGD': 0.0, '24LGF': 0.0, '24LGG': 0.0, '24LGCRF01_AE': 0.0, '24LGCRF02_AG': 0.0, '24LGOther': 0.0, '24LGAll': 0.0, '25DNA': 0.0, '25DNB': 0.3, '25DNC': 0.1, '25DND': 0.0, '25DNF': 0.1, '25DNG': 0.0, '25DNCRF01_AE': 0.3, '25DNCRF02_AG': 0.0, '25DNOther': 0.1, '25DNAll': 0.2, '25DGA': 0.3, '25DGB': 0.0, '25DGC': 0.0, '25DGD': 0.0, '25DGF': 0.0, '25DGG': 0.0, '25DGCRF01_AE': 0.0, '25DGCRF02_AG': 0.0, '25DGOther': 0.1, '25DGAll': 0.0, '25DEA': 0.0, '25DEB': 0.0, '25DEC': 0.1, '25DED': 0.0, '25DEF': 0.0, '25DEG': 0.0, '25DECRF01_AE': 0.0, '25DECRF02_AG': 0.0, '25DEOther': 0.0, '25DEAll': 0.0, '25DVA': 0.0, '25DVB': 0.0, '25DVC': 0.0, '25DVD': 0.0, '25DVF': 0.1, '25DVG': 0.0, '25DVCRF01_AE': 0.0, '25DVCRF02_AG': 0.0, '25DVOther': 0.1, '25DVAll': 0.0, '25DYA': 0.0, '25DYB': 0.0, '25DYC': 0.0, '25DYD': 0.0, '25DYF': 0.1, '25DYG': 0.0, '25DYCRF01_AE': 0.0, '25DYCRF02_AG': 0.0, '25DYOther': 0.0, '25DYAll': 0.0, '25DHA': 0.0, '25DHB': 0.0, '25DHC': 0.0, '25DHD': 0.0, '25DHF': 0.0, '25DHG': 0.0, '25DHCRF01_AE': 0.0, '25DHCRF02_AG': 0.0, '25DHOther': 0.0, '25DHAll': 0.0, '25DAA': 0.0, '25DAB': 0.0, '25DAC': 0.0, '25DAD': 0.0, '25DAF': 0.0, '25DAG': 0.0, '25DACRF01_AE': 0.0, '25DACRF02_AG': 0.0, '25DAOther': 0.0, '25DAAll': 0.0, '25DFA': 0.0, '25DFB': 0.0, '25DFC': 0.0, '25DFD': 0.0, '25DFF': 0.0, '25DFG': 0.0, '25DFCRF01_AE': 0.0, '25DFCRF02_AG': 0.0, '25DFOther': 0.0, '25DFAll': 0.0, '25DPA': 0.0, '25DPB': 0.0, '25DPC': 0.0, '25DPD': 0.0, '25DPF': 0.0, '25DPG': 0.0, '25DPCRF01_AE': 0.0, '25DPCRF02_AG': 0.0, '25DPOther': 0.0, '25DPAll': 0.0, '26TPA': 0.0, '26TPB': 0.0, '26TPC': 0.0, '26TPD': 0.0, '26TPF': 0.1, '26TPG': 0.0, '26TPCRF01_AE': 0.0, '26TPCRF02_AG': 0.0, '26TPOther': 0.3, '26TPAll': 0.1, '26TSA': 0.0, '26TSB': 0.0, '26TSC': 0.0, '26TSD': 0.0, '26TSF': 0.0, '26TSG': 0.0, '26TSCRF01_AE': 0.1, '26TSCRF02_AG': 0.0, '26TSOther': 0.0, '26TSAll': 0.0, '26TAA': 0.0, '26TAB': 0.0, '26TAC': 0.0, '26TAD': 0.0, '26TAF': 0.0, '26TAG': 0.0, '26TACRF01_AE': 0.0, '26TACRF02_AG': 0.0, '26TAOther': 0.1, '26TAAll': 0.0, '26TRA': 0.0, '26TRB': 0.0, '26TRC': 0.0, '26TRD': 0.0, '26TRF': 0.0, '26TRG': 0.0, '26TRCRF01_AE': 0.0, '26TRCRF02_AG': 0.0, '26TROther': 0.1, '26TRAll': 0.0, '26TKA': 0.0, '26TKB': 0.0, '26TKC': 0.0, '26TKD': 0.0, '26TKF': 0.0, '26TKG': 0.0, '26TKCRF01_AE': 0.0, '26TKCRF02_AG': 0.0, '26TKOther': 0.0, '26TKAll': 0.0, '26TQA': 0.0, '26TQB': 0.0, '26TQC': 0.0, '26TQD': 0.0, '26TQF': 0.0, '26TQG': 0.0, '26TQCRF01_AE': 0.0, '26TQCRF02_AG': 0.0, '26TQOther': 0.0, '26TQAll': 0.0, '26TDA': 0.0, '26TDB': 0.0, '26TDC': 0.0, '26TDD': 0.0, '26TDF': 0.0, '26TDG': 0.0, '26TDCRF01_AE': 0.0, '26TDCRF02_AG': 0.0, '26TDOther': 0.0, '26TDAll': 0.0, '26TIA': 0.0, '26TIB': 0.0, '26TIC': 0.0, '26TID': 0.0, '26TIF': 0.0, '26TIG': 0.0, '26TICRF01_AE': 0.0, '26TICRF02_AG': 0.0, '26TIOther': 0.0, '26TIAll': 0.0, '27GRA': 0.0, '27GRB': 0.0, '27GRC': 0.1, '27GRD': 0.0, '27GRF': 0.0, '27GRG': 0.0, '27GRCRF01_AE': 0.0, '27GRCRF02_AG': 0.0, '27GROther': 0.0, '27GRAll': 0.0, '27GEA': 0.0, '27GEB': 0.0, '27GEC': 0.0, '27GED': 0.0, '27GEF': 0.0, '27GEG': 0.0, '27GECRF01_AE': 0.0, '27GECRF02_AG': 0.0, '27GEOther': 0.1, '27GEAll': 0.0, '27GPA': 0.0, '27GPB': 0.0, '27GPC': 0.0, '27GPD': 0.0, '27GPF': 0.0, '27GPG': 0.0, '27GPCRF01_AE': 0.0, '27GPCRF02_AG': 0.0, '27GPOther': 0.1, '27GPAll': 0.0, '27GVA': 0.0, '27GVB': 0.0, '27GVC': 0.0, '27GVD': 0.0, '27GVF': 0.0, '27GVG': 0.0, '27GVCRF01_AE': 0.0, '27GVCRF02_AG': 0.0, '27GVOther': 0.0, '27GVAll': 0.0, '27GAA': 0.0, '27GAB': 0.0, '27GAC': 0.0, '27GAD': 0.0, '27GAF': 0.0, '27GAG': 0.0, '27GACRF01_AE': 0.0, '27GACRF02_AG': 0.0, '27GAOther': 0.0, '27GAAll': 0.0, '27GWA': 0.0, '27GWB': 0.0, '27GWC': 0.0, '27GWD': 0.0, '27GWF': 0.0, '27GWG': 0.0, '27GWCRF01_AE': 0.0, '27GWCRF02_AG': 0.0, '27GWOther': 0.0, '27GWAll': 0.0, '28AQA': 0.0, '28AQB': 0.0, '28AQC': 0.0, '28AQD': 0.0, '28AQF': 0.0, '28AQG': 0.3, '28AQCRF01_AE': 0.0, '28AQCRF02_AG': 0.0, '28AQOther': 0.0, '28AQAll': 0.0, '28APA': 0.0, '28APB': 0.0, '28APC': 0.0, '28APD': 0.0, '28APF': 0.0, '28APG': 0.0, '28APCRF01_AE': 0.0, '28APCRF02_AG': 0.0, '28APOther': 0.1, '28APAll': 0.0, '28ASA': 0.0, '28ASB': 0.0, '28ASC': 0.0, '28ASD': 0.0, '28ASF': 0.0, '28ASG': 0.0, '28ASCRF01_AE': 0.0, '28ASCRF02_AG': 0.0, '28ASOther': 0.0, '28ASAll': 0.0, '28ATA': 0.0, '28ATB': 0.0, '28ATC': 0.0, '28ATD': 0.0, '28ATF': 0.0, '28ATG': 0.0, '28ATCRF01_AE': 0.0, '28ATCRF02_AG': 0.0, '28ATOther': 0.0, '28ATAll': 0.0, '28AEA': 0.0, '28AEB': 0.0, '28AEC': 0.0, '28AED': 0.0, '28AEF': 0.0, '28AEG': 0.0, '28AECRF01_AE': 0.0, '28AECRF02_AG': 0.0, '28AEOther': 0.0, '28AEAll': 0.0, '28AVA': 0.0, '28AVB': 0.0, '28AVC': 0.0, '28AVD': 0.0, '28AVF': 0.0, '28AVG': 0.0, '28AVCRF01_AE': 0.0, '28AVCRF02_AG': 0.0, '28AVOther': 0.0, '28AVAll': 0.0, '28AGA': 0.0, '28AGB': 0.0, '28AGC': 0.0, '28AGD': 0.0, '28AGF': 0.0, '28AGG': 0.0, '28AGCRF01_AE': 0.0, '28AGCRF02_AG': 0.0, '28AGOther': 0.0, '28AGAll': 0.0, '29DNA': 0.0, '29DNB': 0.0, '29DNC': 0.0, '29DND': 0.0, '29DNF': 0.0, '29DNG': 0.0, '29DNCRF01_AE': 0.0, '29DNCRF02_AG': 0.0, '29DNOther': 0.0, '29DNAll': 0.0, '29DVA': 0.0, '29DVB': 0.1, '29DVC': 0.1, '29DVD': 0.0, '29DVF': 0.1, '29DVG': 0.0, '29DVCRF01_AE': 0.0, '29DVCRF02_AG': 0.0, '29DVOther': 0.0, '29DVAll': 0.1, '29DGA': 0.0, '29DGB': 0.0, '29DGC': 0.1, '29DGD': 0.0, '29DGF': 0.0, '29DGG': 0.0, '29DGCRF01_AE': 0.0, '29DGCRF02_AG': 0.0, '29DGOther': 0.0, '29DGAll': 0.0, '29DSA': 0.0, '29DSB': 0.0, '29DSC': 0.0, '29DSD': 0.0, '29DSF': 0.0, '29DSG': 0.0, '29DSCRF01_AE': 0.0, '29DSCRF02_AG': 0.0, '29DSOther': 0.0, '29DSAll': 0.0, '29DYA': 0.0, '29DYB': 0.0, '29DYC': 0.0, '29DYD': 0.0, '29DYF': 0.0, '29DYG': 0.0, '29DYCRF01_AE': 0.0, '29DYCRF02_AG': 0.0, '29DYOther': 0.0, '29DYAll': 0.0, '29DEA': 0.0, '29DEB': 0.0, '29DEC': 0.0, '29DED': 0.0, '29DEF': 0.0, '29DEG': 0.0, '29DECRF01_AE': 0.0, '29DECRF02_AG': 0.0, '29DEOther': 0.0, '29DEAll': 0.0, '29DAA': 0.0, '29DAB': 0.0, '29DAC': 0.0, '29DAD': 0.0, '29DAF': 0.0, '29DAG': 0.0, '29DACRF01_AE': 0.0, '29DACRF02_AG': 0.0, '29DAOther': 0.0, '29DAAll': 0.0, '29DHA': 0.0, '29DHB': 0.0, '29DHC': 0.0, '29DHD': 0.0, '29DHF': 0.0, '29DHG': 0.0, '29DHCRF01_AE': 0.0, '29DHCRF02_AG': 0.0, '29DHOther': 0.0, '29DHAll': 0.0, '30DNA': 0.0, '30DNB': 6.1, '30DNC': 2.2, '30DND': 6.7, '30DNF': 4.2, '30DNG': 0.8, '30DNCRF01_AE': 0.4, '30DNCRF02_AG': 0.0, '30DNOther': 6.1, '30DNAll': 5.2, '30DYA': 0.0, '30DYB': 0.0, '30DYC': 0.0, '30DYD': 0.0, '30DYF': 0.1, '30DYG': 0.0, '30DYCRF01_AE': 0.1, '30DYCRF02_AG': 0.0, '30DYOther': 0.1, '30DYAll': 0.0, '30DVA': 0.3, '30DVB': 0.0, '30DVC': 0.0, '30DVD': 0.0, '30DVF': 0.0, '30DVG': 0.0, '30DVCRF01_AE': 0.0, '30DVCRF02_AG': 0.0, '30DVOther': 0.0, '30DVAll': 0.0, '30DEA': 0.0, '30DEB': 0.0, '30DEC': 0.0, '30DED': 0.0, '30DEF': 0.0, '30DEG': 0.0, '30DECRF01_AE': 0.0, '30DECRF02_AG': 0.0, '30DEOther': 0.0, '30DEAll': 0.0, '30DQA': 0.0, '30DQB': 0.0, '30DQC': 0.0, '30DQD': 0.0, '30DQF': 0.0, '30DQG': 0.0, '30DQCRF01_AE': 0.0, '30DQCRF02_AG': 0.0, '30DQOther': 0.0, '30DQAll': 0.0, '30DAA': 0.0, '30DAB': 0.0, '30DAC': 0.0, '30DAD': 0.0, '30DAF': 0.0, '30DAG': 0.0, '30DACRF01_AE': 0.0, '30DACRF02_AG': 0.0, '30DAOther': 0.0, '30DAAll': 0.0, '30DHA': 0.0, '30DHB': 0.0, '30DHC': 0.0, '30DHD': 0.0, '30DHF': 0.0, '30DHG': 0.0, '30DHCRF01_AE': 0.0, '30DHCRF02_AG': 0.0, '30DHOther': 0.0, '30DHAll': 0.0, '30DGA': 0.0, '30DGB': 0.0, '30DGC': 0.0, '30DGD': 0.0, '30DGF': 0.0, '30DGG': 0.0, '30DGCRF01_AE': 0.0, '30DGCRF02_AG': 0.0, '30DGOther': 0.0, '30DGAll': 0.0, '31TPA': 0.3, '31TPB': 0.1, '31TPC': 0.0, '31TPD': 0.4, '31TPF': 0.1, '31TPG': 0.0, '31TPCRF01_AE': 0.1, '31TPCRF02_AG': 0.0, '31TPOther': 0.1, '31TPAll': 0.1, '31TAA': 0.3, '31TAB': 0.0, '31TAC': 0.1, '31TAD': 0.0, '31TAF': 0.1, '31TAG': 0.0, '31TACRF01_AE': 0.0, '31TACRF02_AG': 0.0, '31TAOther': 0.1, '31TAAll': 0.1, '31TSA': 0.0, '31TSB': 0.0, '31TSC': 0.0, '31TSD': 0.0, '31TSF': 0.1, '31TSG': 0.0, '31TSCRF01_AE': 0.0, '31TSCRF02_AG': 0.0, '31TSOther': 0.1, '31TSAll': 0.0, '31TIA': 0.0, '31TIB': 0.0, '31TIC': 0.0, '31TID': 0.0, '31TIF': 0.0, '31TIG': 0.0, '31TICRF01_AE': 0.1, '31TICRF02_AG': 0.0, '31TIOther': 0.0, '31TIAll': 0.0, '31TNA': 0.0, '31TNB': 0.0, '31TNC': 0.0, '31TND': 0.0, '31TNF': 0.0, '31TNG': 0.0, '31TNCRF01_AE': 0.0, '31TNCRF02_AG': 0.0, '31TNOther': 0.0, '31TNAll': 0.0, '31TKA': 0.0, '31TKB': 0.0, '31TKC': 0.0, '31TKD': 0.0, '31TKF': 0.0, '31TKG': 0.0, '31TKCRF01_AE': 0.0, '31TKCRF02_AG': 0.0, '31TKOther': 0.0, '31TKAll': 0.0, '31T#A': 0.0, '31T#B': 0.0, '31T#C': 0.0, '31T#D': 0.0, '31T#F': 0.0, '31T#G': 0.0, '31T#CRF01_AE': 0.0, '31T#CRF02_AG': 0.0, '31T#Other': 0.0, '31T#All': 0.0, '31TRA': 0.0, '31TRB': 0.0, '31TRC': 0.0, '31TRD': 0.0, '31TRF': 0.0, '31TRG': 0.0, '31TRCRF01_AE': 0.0, '31TRCRF02_AG': 0.0, '31TROther': 0.0, '31TRAll': 0.0, '31TGA': 0.0, '31TGB': 0.0, '31TGC': 0.0, '31TGD': 0.0, '31TGF': 0.0, '31TGG': 0.0, '31TGCRF01_AE': 0.0, '31TGCRF02_AG': 0.0, '31TGOther': 0.0, '31TGAll': 0.0, '32VIA': 1.2, '32VIB': 5.7, '32VIC': 0.6, '32VID': 5.0, '32VIF': 3.0, '32VIG': 1.3, '32VICRF01_AE': 0.7, '32VICRF02_AG': 0.0, '32VIOther': 1.8, '32VIAll': 4.4, '32VAA': 0.0, '32VAB': 0.1, '32VAC': 0.1, '32VAD': 0.4, '32VAF': 0.4, '32VAG': 0.0, '32VACRF01_AE': 0.0, '32VACRF02_AG': 0.0, '32VAOther': 0.1, '32VAAll': 0.1, '32VLA': 0.3, '32VLB': 0.0, '32VLC': 0.0, '32VLD': 0.0, '32VLF': 0.0, '32VLG': 0.0, '32VLCRF01_AE': 0.0, '32VLCRF02_AG': 0.0, '32VLOther': 0.1, '32VLAll': 0.0, '32VEA': 0.0, '32VEB': 0.0, '32VEC': 0.0, '32VED': 0.0, '32VEF': 0.0, '32VEG': 0.0, '32VECRF01_AE': 0.0, '32VECRF02_AG': 0.0, '32VEOther': 0.1, '32VEAll': 0.0, '32VTA': 0.0, '32VTB': 0.0, '32VTC': 0.0, '32VTD': 0.0, '32VTF': 0.0, '32VTG': 0.0, '32VTCRF01_AE': 0.0, '32VTCRF02_AG': 0.0, '32VTOther': 0.0, '32VTAll': 0.0, '32VMA': 0.0, '32VMB': 0.0, '32VMC': 0.0, '32VMD': 0.0, '32VMF': 0.0, '32VMG': 0.0, '32VMCRF01_AE': 0.0, '32VMCRF02_AG': 0.0, '32VMOther': 0.0, '32VMAll': 0.0, '32V#A': 0.0, '32V#B': 0.0, '32V#C': 0.0, '32V#D': 0.0, '32V#F': 0.0, '32V#G': 0.0, '32V#CRF01_AE': 0.0, '32V#CRF02_AG': 0.0, '32V#Other': 0.0, '32V#All': 0.0, '32VPA': 0.0, '32VPB': 0.0, '32VPC': 0.0, '32VPD': 0.0, '32VPF': 0.0, '32VPG': 0.0, '32VPCRF01_AE': 0.0, '32VPCRF02_AG': 0.0, '32VPOther': 0.0, '32VPAll': 0.0, '32VHA': 0.0, '32VHB': 0.0, '32VHC': 0.0, '32VHD': 0.0, '32VHF': 0.0, '32VHG': 0.0, '32VHCRF01_AE': 0.0, '32VHCRF02_AG': 0.0, '32VHOther': 0.0, '32VHAll': 0.0, '32VGA': 0.0, '32VGB': 0.0, '32VGC': 0.0, '32VGD': 0.0, '32VGF': 0.0, '32VGG': 0.0, '32VGCRF01_AE': 0.0, '32VGCRF02_AG': 0.0, '32VGOther': 0.0, '32VGAll': 0.0, '33LFA': 6.8, '33LFB': 13.5, '33LFC': 4.3, '33LFD': 12.2, '33LFF': 12.9, '33LFG': 2.3, '33LFCRF01_AE': 4.6, '33LFCRF02_AG': 1.5, '33LFOther': 7.5, '33LFAll': 11.3, '33LIA': 3.9, '33LIB': 2.4, '33LIC': 0.4, '33LID': 1.7, '33LIF': 2.5, '33LIG': 2.9, '33LICRF01_AE': 2.7, '33LICRF02_AG': 4.2, '33LIOther': 0.9, '33LIAll': 2.1, '33LVA': 0.3, '33LVB': 2.1, '33LVC': 0.4, '33LVD': 2.9, '33LVF': 0.3, '33LVG': 0.8, '33LVCRF01_AE': 0.6, '33LVCRF02_AG': 1.1, '33LVOther': 0.0, '33LVAll': 1.5, '33L#A': 0.3, '33L#B': 0.2, '33L#C': 0.1, '33L#D': 0.4, '33L#F': 0.1, '33L#G': 0.3, '33L#CRF01_AE': 0.0, '33L#CRF02_AG': 0.4, '33L#Other': 0.2, '33L#All': 0.2, '33LMA': 0.0, '33LMB': 0.2, '33LMC': 0.0, '33LMD': 0.0, '33LMF': 0.0, '33LMG': 0.0, '33LMCRF01_AE': 0.1, '33LMCRF02_AG': 0.0, '33LMOther': 0.0, '33LMAll': 0.1, '33LSA': 0.0, '33LSB': 0.0, '33LSC': 0.0, '33LSD': 0.0, '33LSF': 0.0, '33LSG': 0.0, '33LSCRF01_AE': 0.0, '33LSCRF02_AG': 0.0, '33LSOther': 0.0, '33LSAll': 0.0, '33LTA': 0.0, '33LTB': 0.0, '33LTC': 0.0, '33LTD': 0.0, '33LTF': 0.0, '33LTG': 0.0, '33LTCRF01_AE': 0.0, '33LTCRF02_AG': 0.0, '33LTOther': 0.0, '33LTAll': 0.0, '33LKA': 0.0, '33LKB': 0.0, '33LKC': 0.0, '33LKD': 0.0, '33LKF': 0.0, '33LKG': 0.0, '33LKCRF01_AE': 0.0, '33LKCRF02_AG': 0.0, '33LKOther': 0.0, '33LKAll': 0.0, '34EQA': 0.0, '34EQB': 3.4, '34EQC': 0.4, '34EQD': 0.4, '34EQF': 0.5, '34EQG': 0.5, '34EQCRF01_AE': 0.1, '34EQCRF02_AG': 0.0, '34EQOther': 0.5, '34EQAll': 2.4, '34EKA': 0.3, '34EKB': 0.5, '34EKC': 0.2, '34EKD': 0.8, '34EKF': 0.4, '34EKG': 0.0, '34EKCRF01_AE': 0.0, '34EKCRF02_AG': 0.0, '34EKOther': 0.5, '34EKAll': 0.4, '34EAA': 0.0, '34EAB': 0.1, '34EAC': 0.0, '34EAD': 0.4, '34EAF': 0.5, '34EAG': 0.0, '34EACRF01_AE': 0.1, '34EACRF02_AG': 0.0, '34EAOther': 0.3, '34EAAll': 0.1, '34EDA': 0.0, '34EDB': 0.4, '34EDC': 0.1, '34EDD': 0.0, '34EDF': 0.1, '34EDG': 0.3, '34EDCRF01_AE': 0.0, '34EDCRF02_AG': 0.0, '34EDOther': 0.2, '34EDAll': 0.3, '34EGA': 0.3, '34EGB': 0.0, '34EGC': 0.0, '34EGD': 0.0, '34EGF': 0.1, '34EGG': 0.0, '34EGCRF01_AE': 0.1, '34EGCRF02_AG': 0.4, '34EGOther': 0.0, '34EGAll': 0.0, '34ENA': 0.3, '34ENB': 0.2, '34ENC': 0.0, '34END': 0.4, '34ENF': 0.5, '34ENG': 0.0, '34ENCRF01_AE': 0.0, '34ENCRF02_AG': 0.0, '34ENOther': 0.1, '34ENAll': 0.2, '34ETA': 0.0, '34ETB': 0.1, '34ETC': 0.0, '34ETD': 0.0, '34ETF': 0.0, '34ETG': 0.0, '34ETCRF01_AE': 0.0, '34ETCRF02_AG': 0.0, '34ETOther': 0.1, '34ETAll': 0.1, '34EVA': 0.0, '34EVB': 0.1, '34EVC': 0.0, '34EVD': 0.0, '34EVF': 0.1, '34EVG': 0.0, '34EVCRF01_AE': 0.1, '34EVCRF02_AG': 0.0, '34EVOther': 0.2, '34EVAll': 0.1, '34ERA': 0.0, '34ERB': 0.1, '34ERC': 0.0, '34ERD': 0.0, '34ERF': 0.1, '34ERG': 0.0, '34ERCRF01_AE': 0.0, '34ERCRF02_AG': 0.0, '34EROther': 0.1, '34ERAll': 0.1, '34ESA': 0.0, '34ESB': 0.0, '34ESC': 0.0, '34ESD': 0.0, '34ESF': 0.0, '34ESG': 0.0, '34ESCRF01_AE': 0.0, '34ESCRF02_AG': 0.0, '34ESOther': 0.0, '34ESAll': 0.0, '34ELA': 0.0, '34ELB': 0.0, '34ELC': 0.0, '34ELD': 0.0, '34ELF': 0.0, '34ELG': 0.0, '34ELCRF01_AE': 0.0, '34ELCRF02_AG': 0.0, '34ELOther': 0.0, '34ELAll': 0.0, '35EDA': 90.2, '35EDB': 33.0, '35EDC': 30.6, '35EDD': 50.6, '35EDF': 89.6, '35EDG': 71.3, '35EDCRF01_AE': 88.0, '35EDCRF02_AG': 26.9, '35EDOther': 83.4, '35EDAll': 43.0, '35ENA': 3.6, '35ENB': 0.8, '35ENC': 0.7, '35END': 1.7, '35ENF': 4.3, '35ENG': 8.4, '35ENCRF01_AE': 2.1, '35ENCRF02_AG': 3.8, '35ENOther': 2.9, '35ENAll': 1.4, '35EQA': 0.0, '35EQB': 0.1, '35EQC': 0.1, '35EQD': 0.0, '35EQF': 0.0, '35EQG': 5.7, '35EQCRF01_AE': 0.0, '35EQCRF02_AG': 1.5, '35EQOther': 0.1, '35EQAll': 0.2, '35EGA': 0.6, '35EGB': 0.9, '35EGC': 0.9, '35EGD': 2.1, '35EGF': 0.2, '35EGG': 2.3, '35EGCRF01_AE': 0.3, '35EGCRF02_AG': 5.8, '35EGOther': 1.0, '35EGAll': 0.9, '35EKA': 0.0, '35EKB': 0.1, '35EKC': 0.0, '35EKD': 0.4, '35EKF': 0.1, '35EKG': 0.3, '35EKCRF01_AE': 0.0, '35EKCRF02_AG': 1.9, '35EKOther': 0.0, '35EKAll': 0.1, '35EHA': 0.3, '35EHB': 0.0, '35EHC': 0.0, '35EHD': 0.4, '35EHF': 0.1, '35EHG': 1.0, '35EHCRF01_AE': 0.0, '35EHCRF02_AG': 0.0, '35EHOther': 0.1, '35EHAll': 0.0, '35E#A': 0.3, '35E#B': 0.1, '35E#C': 0.3, '35E#D': 0.0, '35E#F': 0.1, '35E#G': 0.0, '35E#CRF01_AE': 0.0, '35E#CRF02_AG': 0.0, '35E#Other': 0.1, '35E#All': 0.1, '35EAA': 0.0, '35EAB': 0.0, '35EAC': 0.0, '35EAD': 0.0, '35EAF': 0.0, '35EAG': 0.0, '35EACRF01_AE': 0.0, '35EACRF02_AG': 0.8, '35EAOther': 0.0, '35EAAll': 0.0, '35ESA': 0.0, '35ESB': 0.0, '35ESC': 0.0, '35ESD': 0.0, '35ESF': 0.1, '35ESG': 0.8, '35ESCRF01_AE': 0.1, '35ESCRF02_AG': 0.0, '35ESOther': 0.0, '35ESAll': 0.0, '35ERA': 0.0, '35ERB': 0.0, '35ERC': 0.0, '35ERD': 0.0, '35ERF': 0.0, '35ERG': 0.3, '35ERCRF01_AE': 0.0, '35ERCRF02_AG': 0.0, '35EROther': 0.0, '35ERAll': 0.0, '35ETA': 0.0, '35ETB': 0.0, '35ETC': 0.0, '35ETD': 0.0, '35ETF': 0.0, '35ETG': 0.0, '35ETCRF01_AE': 0.0, '35ETCRF02_AG': 0.0, '35ETOther': 0.0, '35ETAll': 0.0, '35EYA': 0.0, '35EYB': 0.0, '35EYC': 0.0, '35EYD': 0.0, '35EYF': 0.0, '35EYG': 0.0, '35EYCRF01_AE': 0.0, '35EYCRF02_AG': 0.0, '35EYOther': 0.0, '35EYAll': 0.0, '35EVA': 0.0, '35EVB': 0.0, '35EVC': 0.0, '35EVD': 0.0, '35EVF': 0.0, '35EVG': 0.0, '35EVCRF01_AE': 0.0, '35EVCRF02_AG': 0.0, '35EVOther': 0.0, '35EVAll': 0.0, '35EIA': 0.0, '35EIB': 0.0, '35EIC': 0.0, '35EID': 0.0, '35EIF': 0.0, '35EIG': 0.0, '35EICRF01_AE': 0.0, '35EICRF02_AG': 0.0, '35EIOther': 0.0, '35EIAll': 0.0, '36MIA': 97.6, '36MIB': 35.7, '36MIC': 84.1, '36MID': 78.7, '36MIF': 94.8, '36MIG': 98.2, '36MICRF01_AE': 97.9, '36MICRF02_AG': 99.2, '36MIOther': 95.6, '36MIAll': 54.1, '36MLA': 1.2, '36MLB': 2.9, '36MLC': 3.8, '36MLD': 1.3, '36MLF': 1.4, '36MLG': 0.0, '36MLCRF01_AE': 0.4, '36MLCRF02_AG': 0.0, '36MLOther': 1.2, '36MLAll': 2.6, '36MVA': 0.6, '36MVB': 1.7, '36MVC': 1.1, '36MVD': 0.4, '36MVF': 2.1, '36MVG': 1.3, '36MVCRF01_AE': 0.4, '36MVCRF02_AG': 0.0, '36MVOther': 1.6, '36MVAll': 1.4, '36MTA': 0.0, '36MTB': 0.2, '36MTC': 1.0, '36MTD': 0.8, '36MTF': 0.1, '36MTG': 0.0, '36MTCRF01_AE': 0.0, '36MTCRF02_AG': 0.4, '36MTOther': 0.2, '36MTAll': 0.3, '36M#A': 0.0, '36M#B': 0.1, '36M#C': 0.0, '36M#D': 0.0, '36M#F': 0.0, '36M#G': 0.3, '36M#CRF01_AE': 0.0, '36M#CRF02_AG': 0.0, '36M#Other': 0.0, '36M#All': 0.0, '36MAA': 0.0, '36MAB': 0.1, '36MAC': 0.1, '36MAD': 0.0, '36MAF': 0.0, '36MAG': 0.0, '36MACRF01_AE': 0.0, '36MACRF02_AG': 0.0, '36MAOther': 0.0, '36MAAll': 0.1, '36MDA': 0.0, '36MDB': 0.0, '36MDC': 0.0, '36MDD': 0.0, '36MDF': 0.1, '36MDG': 0.0, '36MDCRF01_AE': 0.0, '36MDCRF02_AG': 0.0, '36MDOther': 0.1, '36MDAll': 0.0, '36MKA': 0.0, '36MKB': 0.0, '36MKC': 0.0, '36MKD': 0.0, '36MKF': 0.0, '36MKG': 0.0, '36MKCRF01_AE': 0.0, '36MKCRF02_AG': 0.0, '36MKOther': 0.0, '36MKAll': 0.0, '36MSA': 0.0, '36MSB': 0.0, '36MSC': 0.0, '36MSD': 0.0, '36MSF': 0.0, '36MSG': 0.0, '36MSCRF01_AE': 0.0, '36MSCRF02_AG': 0.0, '36MSOther': 0.0, '36MSAll': 0.0, '36MFA': 0.0, '36MFB': 0.0, '36MFC': 0.0, '36MFD': 0.0, '36MFF': 0.0, '36MFG': 0.0, '36MFCRF01_AE': 0.0, '36MFCRF02_AG': 0.0, '36MFOther': 0.0, '36MFAll': 0.0, '36MNA': 0.0, '36MNB': 0.0, '36MNC': 0.0, '36MND': 0.0, '36MNF': 0.0, '36MNG': 0.0, '36MNCRF01_AE': 0.0, '36MNCRF02_AG': 0.0, '36MNOther': 0.0, '36MNAll': 0.0, '36MEA': 0.0, '36MEB': 0.0, '36MEC': 0.0, '36MED': 0.0, '36MEF': 0.0, '36MEG': 0.0, '36MECRF01_AE': 0.0, '36MECRF02_AG': 0.0, '36MEOther': 0.0, '36MEAll': 0.0, '36MQA': 0.0, '36MQB': 0.0, '36MQC': 0.0, '36MQD': 0.0, '36MQF': 0.0, '36MQG': 0.0, '36MQCRF01_AE': 0.0, '36MQCRF02_AG': 0.0, '36MQOther': 0.0, '36MQAll': 0.0, '36MHA': 0.0, '36MHB': 0.0, '36MHC': 0.0, '36MHD': 0.0, '36MHF': 0.0, '36MHG': 0.0, '36MHCRF01_AE': 0.0, '36MHCRF02_AG': 0.0, '36MHOther': 0.0, '36MHAll': 0.0, '36MRA': 0.0, '36MRB': 0.0, '36MRC': 0.0, '36MRD': 0.0, '36MRF': 0.0, '36MRG': 0.0, '36MRCRF01_AE': 0.0, '36MRCRF02_AG': 0.0, '36MROther': 0.0, '36MRAll': 0.0, '36MGA': 0.0, '36MGB': 0.0, '36MGC': 0.0, '36MGD': 0.0, '36MGF': 0.0, '36MGG': 0.0, '36MGCRF01_AE': 0.0, '36MGCRF02_AG': 0.0, '36MGOther': 0.0, '36MGAll': 0.0, '37NDA': 19.9, '37NDB': 15.7, '37NDC': 6.8, '37NDD': 13.0, '37NDF': 14.3, '37NDG': 18.0, '37NDCRF01_AE': 8.3, '37NDCRF02_AG': 21.8, '37NDOther': 17.7, '37NDAll': 14.7, '37NSA': 0.6, '37NSB': 12.0, '37NSC': 6.4, '37NSD': 5.5, '37NSF': 1.2, '37NSG': 1.3, '37NSCRF01_AE': 0.6, '37NSCRF02_AG': 6.1, '37NSOther': 1.4, '37NSAll': 9.1, '37NEA': 3.0, '37NEB': 6.4, '37NEC': 1.7, '37NED': 6.3, '37NEF': 2.1, '37NEG': 2.9, '37NECRF01_AE': 0.7, '37NECRF02_AG': 3.1, '37NEOther': 3.6, '37NEAll': 5.1, '37NKA': 0.0, '37NKB': 0.3, '37NKC': 29.0, '37NKD': 0.4, '37NKF': 0.9, '37NKG': 0.8, '37NKCRF01_AE': 0.0, '37NKCRF02_AG': 0.4, '37NKOther': 0.6, '37NKAll': 3.7, '37NTA': 0.6, '37NTB': 3.0, '37NTC': 0.8, '37NTD': 2.1, '37NTF': 3.3, '37NTG': 0.5, '37NTCRF01_AE': 0.1, '37NTCRF02_AG': 1.1, '37NTOther': 1.2, '37NTAll': 2.4, '37NHA': 0.3, '37NHB': 1.2, '37NHC': 0.1, '37NHD': 0.4, '37NHF': 0.3, '37NHG': 0.0, '37NHCRF01_AE': 0.3, '37NHCRF02_AG': 0.4, '37NHOther': 0.4, '37NHAll': 0.9, '37NCA': 0.0, '37NCB': 1.2, '37NCC': 0.3, '37NCD': 0.4, '37NCF': 0.0, '37NCG': 0.0, '37NCCRF01_AE': 0.0, '37NCCRF02_AG': 0.0, '37NCOther': 0.0, '37NCAll': 0.8, '37NAA': 0.0, '37NAB': 0.8, '37NAC': 0.3, '37NAD': 0.0, '37NAF': 0.0, '37NAG': 0.0, '37NACRF01_AE': 0.0, '37NACRF02_AG': 0.0, '37NAOther': 0.1, '37NAAll': 0.6, '37NQA': 0.0, '37NQB': 0.6, '37NQC': 0.8, '37NQD': 0.4, '37NQF': 0.1, '37NQG': 0.0, '37NQCRF01_AE': 0.0, '37NQCRF02_AG': 0.0, '37NQOther': 0.0, '37NQAll': 0.5, '37NYA': 0.0, '37NYB': 0.4, '37NYC': 0.0, '37NYD': 0.4, '37NYF': 0.1, '37NYG': 0.0, '37NYCRF01_AE': 0.1, '37NYCRF02_AG': 0.4, '37NYOther': 0.1, '37NYAll': 0.3, '37NIA': 0.0, '37NIB': 0.2, '37NIC': 0.0, '37NID': 0.0, '37NIF': 0.5, '37NIG': 0.3, '37NICRF01_AE': 0.0, '37NICRF02_AG': 0.0, '37NIOther': 0.0, '37NIAll': 0.1, '37NRA': 0.0, '37NRB': 0.1, '37NRC': 1.0, '37NRD': 0.0, '37NRF': 0.0, '37NRG': 0.0, '37NRCRF01_AE': 0.0, '37NRCRF02_AG': 0.0, '37NROther': 0.0, '37NRAll': 0.2, '37N#A': 0.0, '37N#B': 0.0, '37N#C': 0.0, '37N#D': 0.0, '37N#F': 0.1, '37N#G': 0.0, '37N#CRF01_AE': 0.0, '37N#CRF02_AG': 0.4, '37N#Other': 0.1, '37N#All': 0.0, '37NPA': 0.0, '37NPB': 0.3, '37NPC': 0.0, '37NPD': 0.0, '37NPF': 0.0, '37NPG': 0.0, '37NPCRF01_AE': 0.0, '37NPCRF02_AG': 0.0, '37NPOther': 0.0, '37NPAll': 0.2, '37NGA': 0.0, '37NGB': 0.0, '37NGC': 0.0, '37NGD': 0.4, '37NGF': 0.1, '37NGG': 0.0, '37NGCRF01_AE': 0.0, '37NGCRF02_AG': 0.0, '37NGOther': 0.1, '37NGAll': 0.0, '37NFA': 0.0, '37NFB': 0.0, '37NFC': 0.0, '37NFD': 0.0, '37NFF': 0.0, '37NFG': 0.0, '37NFCRF01_AE': 0.0, '37NFCRF02_AG': 0.0, '37NFOther': 0.1, '37NFAll': 0.0, '37NVA': 0.0, '37NVB': 0.0, '37NVC': 0.0, '37NVD': 0.0, '37NVF': 0.0, '37NVG': 0.0, '37NVCRF01_AE': 0.0, '37NVCRF02_AG': 0.0, '37NVOther': 0.1, '37NVAll': 0.0, '37NMA': 0.0, '37NMB': 0.0, '37NMC': 0.0, '37NMD': 0.0, '37NMF': 0.0, '37NMG': 0.0, '37NMCRF01_AE': 0.0, '37NMCRF02_AG': 0.0, '37NMOther': 0.0, '37NMAll': 0.0, '37NWA': 0.0, '37NWB': 0.0, '37NWC': 0.0, '37NWD': 0.0, '37NWF': 0.0, '37NWG': 0.0, '37NWCRF01_AE': 0.0, '37NWCRF02_AG': 0.0, '37NWOther': 0.0, '37NWAll': 0.0, '38LFA': 0.0, '38LFB': 0.2, '38LFC': 0.0, '38LFD': 0.0, '38LFF': 0.1, '38LFG': 0.5, '38LFCRF01_AE': 0.0, '38LFCRF02_AG': 0.4, '38LFOther': 0.1, '38LFAll': 0.1, '38LIA': 0.0, '38LIB': 0.2, '38LIC': 0.1, '38LID': 0.0, '38LIF': 0.0, '38LIG': 0.0, '38LICRF01_AE': 0.0, '38LICRF02_AG': 0.4, '38LIOther': 0.0, '38LIAll': 0.1, '38LMA': 0.0, '38LMB': 0.1, '38LMC': 0.2, '38LMD': 0.0, '38LMF': 0.0, '38LMG': 0.0, '38LMCRF01_AE': 0.0, '38LMCRF02_AG': 0.0, '38LMOther': 0.1, '38LMAll': 0.1, '38LWA': 0.0, '38LWB': 0.2, '38LWC': 0.0, '38LWD': 0.0, '38LWF': 0.1, '38LWG': 0.0, '38LWCRF01_AE': 0.0, '38LWCRF02_AG': 0.0, '38LWOther': 0.2, '38LWAll': 0.2, '38LVA': 0.0, '38LVB': 0.1, '38LVC': 0.1, '38LVD': 0.0, '38LVF': 0.1, '38LVG': 0.0, '38LVCRF01_AE': 0.1, '38LVCRF02_AG': 0.0, '38LVOther': 0.1, '38LVAll': 0.1, '38LKA': 0.0, '38LKB': 0.0, '38LKC': 0.0, '38LKD': 0.0, '38LKF': 0.0, '38LKG': 0.0, '38LKCRF01_AE': 0.1, '38LKCRF02_AG': 0.0, '38LKOther': 0.0, '38LKAll': 0.0, '38LSA': 0.0, '38LSB': 0.0, '38LSC': 0.0, '38LSD': 0.0, '38LSF': 0.0, '38LSG': 0.0, '38LSCRF01_AE': 0.0, '38LSCRF02_AG': 0.0, '38LSOther': 0.1, '38LSAll': 0.0, '38LYA': 0.0, '38LYB': 0.0, '38LYC': 0.0, '38LYD': 0.0, '38LYF': 0.0, '38LYG': 0.0, '38LYCRF01_AE': 0.0, '38LYCRF02_AG': 0.0, '38LYOther': 0.0, '38LYAll': 0.0, '38LCA': 0.0, '38LCB': 0.0, '38LCC': 0.0, '38LCD': 0.0, '38LCF': 0.0, '38LCG': 0.0, '38LCCRF01_AE': 0.0, '38LCCRF02_AG': 0.0, '38LCOther': 0.0, '38LCAll': 0.0, '38L#A': 0.0, '38L#B': 0.0, '38L#C': 0.0, '38L#D': 0.0, '38L#F': 0.0, '38L#G': 0.0, '38L#CRF01_AE': 0.0, '38L#CRF02_AG': 0.0, '38L#Other': 0.0, '38L#All': 0.0, '38LPA': 0.0, '38LPB': 0.0, '38LPC': 0.0, '38LPD': 0.0, '38LPF': 0.0, '38LPG': 0.0, '38LPCRF01_AE': 0.0, '38LPCRF02_AG': 0.0, '38LPOther': 0.0, '38LPAll': 0.0, '38LRA': 0.0, '38LRB': 0.0, '38LRC': 0.0, '38LRD': 0.0, '38LRF': 0.0, '38LRG': 0.0, '38LRCRF01_AE': 0.0, '38LRCRF02_AG': 0.0, '38LROther': 0.0, '38LRAll': 0.0, '38LGA': 0.0, '38LGB': 0.0, '38LGC': 0.0, '38LGD': 0.0, '38LGF': 0.0, '38LGG': 0.0, '38LGCRF01_AE': 0.0, '38LGCRF02_AG': 0.0, '38LGOther': 0.0, '38LGAll': 0.0, '39PSA': 2.4, '39PSB': 2.5, '39PSC': 1.7, '39PSD': 2.5, '39PSF': 5.9, '39PSG': 0.8, '39PSCRF01_AE': 1.8, '39PSCRF02_AG': 3.1, '39PSOther': 2.0, '39PSAll': 2.5, '39PQA': 0.9, '39PQB': 0.9, '39PQC': 0.6, '39PQD': 1.7, '39PQF': 3.4, '39PQG': 0.8, '39PQCRF01_AE': 0.3, '39PQCRF02_AG': 1.1, '39PQOther': 1.5, '39PQAll': 1.0, '39PTA': 0.6, '39PTB': 0.4, '39PTC': 0.3, '39PTD': 0.0, '39PTF': 2.4, '39PTG': 0.0, '39PTCRF01_AE': 0.0, '39PTCRF02_AG': 0.4, '39PTOther': 1.0, '39PTAll': 0.6, '39PLA': 0.0, '39PLB': 0.0, '39PLC': 0.1, '39PLD': 0.0, '39PLF': 0.0, '39PLG': 0.0, '39PLCRF01_AE': 0.0, '39PLCRF02_AG': 1.5, '39PLOther': 0.2, '39PLAll': 0.1, '39PAA': 0.0, '39PAB': 0.1, '39PAC': 0.1, '39PAD': 0.0, '39PAF': 0.1, '39PAG': 0.0, '39PACRF01_AE': 0.0, '39PACRF02_AG': 0.4, '39PAOther': 0.3, '39PAAll': 0.1, '39PVA': 0.0, '39PVB': 0.0, '39PVC': 0.0, '39PVD': 0.0, '39PVF': 0.1, '39PVG': 0.0, '39PVCRF01_AE': 0.0, '39PVCRF02_AG': 0.4, '39PVOther': 0.0, '39PVAll': 0.0, '39PEA': 0.0, '39PEB': 0.0, '39PEC': 0.0, '39PED': 0.0, '39PEF': 0.1, '39PEG': 0.0, '39PECRF01_AE': 0.0, '39PECRF02_AG': 0.0, '39PEOther': 0.0, '39PEAll': 0.0, '39PMA': 0.0, '39PMB': 0.0, '39PMC': 0.0, '39PMD': 0.0, '39PMF': 0.0, '39PMG': 0.0, '39PMCRF01_AE': 0.0, '39PMCRF02_AG': 0.0, '39PMOther': 0.1, '39PMAll': 0.0, '39PKA': 0.0, '39PKB': 0.0, '39PKC': 0.0, '39PKD': 0.0, '39PKF': 0.1, '39PKG': 0.0, '39PKCRF01_AE': 0.0, '39PKCRF02_AG': 0.0, '39PKOther': 0.0, '39PKAll': 0.0, '39PIA': 0.0, '39PIB': 0.0, '39PIC': 0.0, '39PID': 0.0, '39PIF': 0.0, '39PIG': 0.0, '39PICRF01_AE': 0.0, '39PICRF02_AG': 0.0, '39PIOther': 0.0, '39PIAll': 0.0, '39P#A': 0.0, '39P#B': 0.0, '39P#C': 0.0, '39P#D': 0.0, '39P#F': 0.0, '39P#G': 0.0, '39P#CRF01_AE': 0.0, '39P#CRF02_AG': 0.0, '39P#Other': 0.0, '39P#All': 0.0, '39PHA': 0.0, '39PHB': 0.0, '39PHC': 0.0, '39PHD': 0.0, '39PHF': 0.0, '39PHG': 0.0, '39PHCRF01_AE': 0.0, '39PHCRF02_AG': 0.0, '39PHOther': 0.0, '39PHAll': 0.0, '39PDA': 0.0, '39PDB': 0.0, '39PDC': 0.0, '39PDD': 0.0, '39PDF': 0.0, '39PDG': 0.0, '39PDCRF01_AE': 0.0, '39PDCRF02_AG': 0.0, '39PDOther': 0.0, '39PDAll': 0.0, '39PRA': 0.0, '39PRB': 0.0, '39PRC': 0.0, '39PRD': 0.0, '39PRF': 0.0, '39PRG': 0.0, '39PRCRF01_AE': 0.0, '39PRCRF02_AG': 0.0, '39PROther': 0.0, '39PRAll': 0.0, '39PGA': 0.0, '39PGB': 0.0, '39PGC': 0.0, '39PGD': 0.0, '39PGF': 0.0, '39PGG': 0.0, '39PGCRF01_AE': 0.0, '39PGCRF02_AG': 0.0, '39PGOther': 0.0, '39PGAll': 0.0, '40GRA': 0.0, '40GRB': 0.0, '40GRC': 0.0, '40GRD': 0.0, '40GRF': 0.1, '40GRG': 0.3, '40GRCRF01_AE': 0.3, '40GRCRF02_AG': 0.0, '40GROther': 0.0, '40GRAll': 0.0, '40GEA': 0.0, '40GEB': 0.0, '40GEC': 0.0, '40GED': 0.0, '40GEF': 0.0, '40GEG': 0.3, '40GECRF01_AE': 0.0, '40GECRF02_AG': 0.0, '40GEOther': 0.0, '40GEAll': 0.0, '40GVA': 0.0, '40GVB': 0.0, '40GVC': 0.0, '40GVD': 0.0, '40GVF': 0.0, '40GVG': 0.0, '40GVCRF01_AE': 0.1, '40GVCRF02_AG': 0.0, '40GVOther': 0.0, '40GVAll': 0.0, '40GTA': 0.0, '40GTB': 0.0, '40GTC': 0.0, '40GTD': 0.0, '40GTF': 0.0, '40GTG': 0.0, '40GTCRF01_AE': 0.0, '40GTCRF02_AG': 0.0, '40GTOther': 0.0, '40GTAll': 0.0, '40GKA': 0.0, '40GKB': 0.0, '40GKC': 0.0, '40GKD': 0.0, '40GKF': 0.0, '40GKG': 0.0, '40GKCRF01_AE': 0.0, '40GKCRF02_AG': 0.0, '40GKOther': 0.0, '40GKAll': 0.0, '40GAA': 0.0, '40GAB': 0.0, '40GAC': 0.0, '40GAD': 0.0, '40GAF': 0.0, '40GAG': 0.0, '40GACRF01_AE': 0.0, '40GACRF02_AG': 0.0, '40GAOther': 0.0, '40GAAll': 0.0, '40G#A': 0.0, '40G#B': 0.0, '40G#C': 0.0, '40G#D': 0.0, '40G#F': 0.0, '40G#G': 0.0, '40G#CRF01_AE': 0.0, '40G#CRF02_AG': 0.0, '40G#Other': 0.0, '40G#All': 0.0, '41RKA': 96.1, '41RKB': 25.8, '41RKC': 61.5, '41RKD': 85.8, '41RKF': 90.7, '41RKG': 92.2, '41RKCRF01_AE': 92.7, '41RKCRF02_AG': 90.4, '41RKOther': 88.1, '41RKAll': 43.6, '41RNA': 0.3, '41RNB': 0.2, '41RNC': 27.8, '41RND': 0.8, '41RNF': 0.3, '41RNG': 0.0, '41RNCRF01_AE': 0.0, '41RNCRF02_AG': 0.4, '41RNOther': 0.1, '41RNAll': 3.4, '41RIA': 0.0, '41RIB': 0.0, '41RIC': 1.5, '41RID': 0.0, '41RIF': 0.1, '41RIG': 0.0, '41RICRF01_AE': 0.0, '41RICRF02_AG': 0.4, '41RIOther': 0.1, '41RIAll': 0.2, '41RTA': 0.0, '41RTB': 0.0, '41RTC': 0.2, '41RTD': 0.8, '41RTF': 0.0, '41RTG': 0.0, '41RTCRF01_AE': 0.0, '41RTCRF02_AG': 0.4, '41RTOther': 0.0, '41RTAll': 0.1, '41RHA': 0.0, '41RHB': 0.0, '41RHC': 0.5, '41RHD': 0.0, '41RHF': 0.0, '41RHG': 0.0, '41RHCRF01_AE': 0.0, '41RHCRF02_AG': 0.0, '41RHOther': 0.0, '41RHAll': 0.1, '41RSA': 0.0, '41RSB': 0.0, '41RSC': 0.2, '41RSD': 0.0, '41RSF': 0.0, '41RSG': 0.0, '41RSCRF01_AE': 0.0, '41RSCRF02_AG': 0.4, '41RSOther': 0.0, '41RSAll': 0.0, '41RPA': 0.0, '41RPB': 0.0, '41RPC': 0.3, '41RPD': 0.0, '41RPF': 0.2, '41RPG': 0.0, '41RPCRF01_AE': 0.0, '41RPCRF02_AG': 0.0, '41RPOther': 0.0, '41RPAll': 0.1, '41REA': 0.0, '41REB': 0.0, '41REC': 0.0, '41RED': 0.0, '41REF': 0.0, '41REG': 0.3, '41RECRF01_AE': 0.0, '41RECRF02_AG': 0.0, '41REOther': 0.2, '41REAll': 0.0, '41RQA': 0.0, '41RQB': 0.0, '41RQC': 0.0, '41RQD': 0.0, '41RQF': 0.0, '41RQG': 0.3, '41RQCRF01_AE': 0.0, '41RQCRF02_AG': 0.0, '41RQOther': 0.1, '41RQAll': 0.0, '41RYA': 0.0, '41RYB': 0.0, '41RYC': 0.1, '41RYD': 0.0, '41RYF': 0.1, '41RYG': 0.0, '41RYCRF01_AE': 0.0, '41RYCRF02_AG': 0.0, '41RYOther': 0.0, '41RYAll': 0.0, '41R#A': 0.0, '41R#B': 0.0, '41R#C': 0.1, '41R#D': 0.0, '41R#F': 0.0, '41R#G': 0.0, '41R#CRF01_AE': 0.0, '41R#CRF02_AG': 0.0, '41R#Other': 0.1, '41R#All': 0.0, '41RFA': 0.0, '41RFB': 0.0, '41RFC': 0.0, '41RFD': 0.0, '41RFF': 0.0, '41RFG': 0.0, '41RFCRF01_AE': 0.0, '41RFCRF02_AG': 0.0, '41RFOther': 0.0, '41RFAll': 0.0, '41RVA': 0.0, '41RVB': 0.0, '41RVC': 0.0, '41RVD': 0.0, '41RVF': 0.0, '41RVG': 0.0, '41RVCRF01_AE': 0.0, '41RVCRF02_AG': 0.0, '41RVOther': 0.0, '41RVAll': 0.0, '41RMA': 0.0, '41RMB': 0.0, '41RMC': 0.0, '41RMD': 0.0, '41RMF': 0.0, '41RMG': 0.0, '41RMCRF01_AE': 0.0, '41RMCRF02_AG': 0.0, '41RMOther': 0.0, '41RMAll': 0.0, '41RLA': 0.0, '41RLB': 0.0, '41RLC': 0.0, '41RLD': 0.0, '41RLF': 0.0, '41RLG': 0.0, '41RLCRF01_AE': 0.0, '41RLCRF02_AG': 0.0, '41RLOther': 0.0, '41RLAll': 0.0, '41RAA': 0.0, '41RAB': 0.0, '41RAC': 0.0, '41RAD': 0.0, '41RAF': 0.0, '41RAG': 0.0, '41RACRF01_AE': 0.0, '41RACRF02_AG': 0.0, '41RAOther': 0.0, '41RAAll': 0.0, '41RDA': 0.0, '41RDB': 0.0, '41RDC': 0.0, '41RDD': 0.0, '41RDF': 0.0, '41RDG': 0.0, '41RDCRF01_AE': 0.0, '41RDCRF02_AG': 0.0, '41RDOther': 0.0, '41RDAll': 0.0, '41RGA': 0.0, '41RGB': 0.0, '41RGC': 0.0, '41RGD': 0.0, '41RGF': 0.0, '41RGG': 0.0, '41RGCRF01_AE': 0.0, '41RGCRF02_AG': 0.0, '41RGOther': 0.0, '41RGAll': 0.0, '42WYA': 0.0, '42WYB': 0.1, '42WYC': 0.0, '42WYD': 0.0, '42WYF': 0.0, '42WYG': 0.0, '42WYCRF01_AE': 0.0, '42WYCRF02_AG': 0.0, '42WYOther': 0.0, '42WYAll': 0.1, '42WGA': 0.3, '42WGB': 0.0, '42WGC': 0.0, '42WGD': 0.0, '42WGF': 0.0, '42WGG': 0.0, '42WGCRF01_AE': 0.0, '42WGCRF02_AG': 0.0, '42WGOther': 0.0, '42WGAll': 0.0, '42WLA': 0.0, '42WLB': 0.0, '42WLC': 0.0, '42WLD': 0.0, '42WLF': 0.0, '42WLG': 0.0, '42WLCRF01_AE': 0.1, '42WLCRF02_AG': 0.0, '42WLOther': 0.0, '42WLAll': 0.0, '42WFA': 0.0, '42WFB': 0.0, '42WFC': 0.0, '42WFD': 0.0, '42WFF': 0.0, '42WFG': 0.0, '42WFCRF01_AE': 0.0, '42WFCRF02_AG': 0.0, '42WFOther': 0.0, '42WFAll': 0.0, '42WSA': 0.0, '42WSB': 0.0, '42WSC': 0.0, '42WSD': 0.0, '42WSF': 0.0, '42WSG': 0.0, '42WSCRF01_AE': 0.0, '42WSCRF02_AG': 0.0, '42WSOther': 0.0, '42WSAll': 0.0, '42WCA': 0.0, '42WCB': 0.0, '42WCC': 0.0, '42WCD': 0.0, '42WCF': 0.0, '42WCG': 0.0, '42WCCRF01_AE': 0.0, '42WCCRF02_AG': 0.0, '42WCOther': 0.0, '42WCAll': 0.0, '42WRA': 0.0, '42WRB': 0.0, '42WRC': 0.0, '42WRD': 0.0, '42WRF': 0.0, '42WRG': 0.0, '42WRCRF01_AE': 0.0, '42WRCRF02_AG': 0.0, '42WROther': 0.0, '42WRAll': 0.0, '43KRA': 2.1, '43KRB': 2.3, '43KRC': 1.3, '43KRD': 15.1, '43KRF': 4.6, '43KRG': 1.8, '43KRCRF01_AE': 1.6, '43KRCRF02_AG': 8.0, '43KROther': 2.6, '43KRAll': 2.5, '43KTA': 1.2, '43KTB': 6.2, '43KTC': 1.8, '43KTD': 5.4, '43KTF': 7.9, '43KTG': 3.1, '43KTCRF01_AE': 0.9, '43KTCRF02_AG': 0.8, '43KTOther': 3.9, '43KTAll': 5.3, '43KIA': 0.3, '43KIB': 0.3, '43KIC': 0.1, '43KID': 0.4, '43KIF': 0.8, '43KIG': 0.3, '43KICRF01_AE': 0.0, '43KICRF02_AG': 0.0, '43KIOther': 0.2, '43KIAll': 0.3, '43KNA': 0.6, '43KNB': 0.2, '43KNC': 0.0, '43KND': 0.4, '43KNF': 0.1, '43KNG': 0.0, '43KNCRF01_AE': 0.0, '43KNCRF02_AG': 0.0, '43KNOther': 0.8, '43KNAll': 0.2, '43KQA': 0.0, '43KQB': 0.2, '43KQC': 0.0, '43KQD': 0.0, '43KQF': 0.3, '43KQG': 0.0, '43KQCRF01_AE': 0.1, '43KQCRF02_AG': 0.0, '43KQOther': 0.2, '43KQAll': 0.2, '43KEA': 0.0, '43KEB': 0.1, '43KEC': 0.1, '43KED': 0.0, '43KEF': 0.2, '43KEG': 0.0, '43KECRF01_AE': 0.0, '43KECRF02_AG': 0.0, '43KEOther': 0.1, '43KEAll': 0.1, '43KSA': 0.0, '43KSB': 0.1, '43KSC': 0.0, '43KSD': 0.0, '43KSF': 0.1, '43KSG': 0.0, '43KSCRF01_AE': 0.0, '43KSCRF02_AG': 0.0, '43KSOther': 0.1, '43KSAll': 0.1, '43KVA': 0.0, '43KVB': 0.0, '43KVC': 0.0, '43KVD': 0.0, '43KVF': 0.0, '43KVG': 0.0, '43KVCRF01_AE': 0.0, '43KVCRF02_AG': 0.0, '43KVOther': 0.0, '43KVAll': 0.0, '43KMA': 0.0, '43KMB': 0.0, '43KMC': 0.0, '43KMD': 0.0, '43KMF': 0.0, '43KMG': 0.0, '43KMCRF01_AE': 0.0, '43KMCRF02_AG': 0.0, '43KMOther': 0.0, '43KMAll': 0.0, '43KPA': 0.0, '43KPB': 0.0, '43KPC': 0.0, '43KPD': 0.0, '43KPF': 0.0, '43KPG': 0.0, '43KPCRF01_AE': 0.0, '43KPCRF02_AG': 0.0, '43KPOther': 0.0, '43KPAll': 0.0, '43KGA': 0.0, '43KGB': 0.0, '43KGC': 0.0, '43KGD': 0.0, '43KGF': 0.0, '43KGG': 0.0, '43KGCRF01_AE': 0.0, '43KGCRF02_AG': 0.0, '43KGOther': 0.0, '43KGAll': 0.0, '44PRA': 0.0, '44PRB': 0.0, '44PRC': 0.0, '44PRD': 0.0, '44PRF': 0.0, '44PRG': 0.0, '44PRCRF01_AE': 0.0, '44PRCRF02_AG': 0.0, '44PROther': 0.0, '44PRAll': 0.0, '44PSA': 0.0, '44PSB': 0.0, '44PSC': 0.0, '44PSD': 0.0, '44PSF': 0.0, '44PSG': 0.0, '44PSCRF01_AE': 0.0, '44PSCRF02_AG': 0.0, '44PSOther': 0.0, '44PSAll': 0.0, '44PTA': 0.0, '44PTB': 0.0, '44PTC': 0.0, '44PTD': 0.0, '44PTF': 0.0, '44PTG': 0.0, '44PTCRF01_AE': 0.0, '44PTCRF02_AG': 0.0, '44PTOther': 0.0, '44PTAll': 0.0, '44PKA': 0.0, '44PKB': 0.0, '44PKC': 0.0, '44PKD': 0.0, '44PKF': 0.0, '44PKG': 0.0, '44PKCRF01_AE': 0.0, '44PKCRF02_AG': 0.0, '44PKOther': 0.0, '44PKAll': 0.0, '44PQA': 0.0, '44PQB': 0.0, '44PQC': 0.0, '44PQD': 0.0, '44PQF': 0.0, '44PQG': 0.0, '44PQCRF01_AE': 0.0, '44PQCRF02_AG': 0.0, '44PQOther': 0.0, '44PQAll': 0.0, '44PLA': 0.0, '44PLB': 0.0, '44PLC': 0.0, '44PLD': 0.0, '44PLF': 0.0, '44PLG': 0.0, '44PLCRF01_AE': 0.0, '44PLCRF02_AG': 0.0, '44PLOther': 0.0, '44PLAll': 0.0, '44PAA': 0.0, '44PAB': 0.0, '44PAC': 0.0, '44PAD': 0.0, '44PAF': 0.0, '44PAG': 0.0, '44PACRF01_AE': 0.0, '44PACRF02_AG': 0.0, '44PAOther': 0.0, '44PAAll': 0.0, '45KRA': 9.8, '45KRB': 3.4, '45KRC': 6.9, '45KRD': 2.9, '45KRF': 7.6, '45KRG': 3.1, '45KRCRF01_AE': 2.7, '45KRCRF02_AG': 1.1, '45KROther': 4.1, '45KRAll': 4.1, '45KVA': 0.6, '45KVB': 0.1, '45KVC': 0.1, '45KVD': 0.8, '45KVF': 0.1, '45KVG': 1.0, '45KVCRF01_AE': 0.3, '45KVCRF02_AG': 0.0, '45KVOther': 0.2, '45KVAll': 0.1, '45KQA': 0.0, '45KQB': 0.3, '45KQC': 0.2, '45KQD': 0.0, '45KQF': 0.5, '45KQG': 0.5, '45KQCRF01_AE': 0.1, '45KQCRF02_AG': 0.0, '45KQOther': 0.5, '45KQAll': 0.3, '45KIA': 0.6, '45KIB': 0.2, '45KIC': 0.1, '45KID': 0.0, '45KIF': 0.3, '45KIG': 0.0, '45KICRF01_AE': 0.1, '45KICRF02_AG': 0.0, '45KIOther': 0.3, '45KIAll': 0.2, '45KTA': 0.0, '45KTB': 0.0, '45KTC': 0.1, '45KTD': 0.8, '45KTF': 0.0, '45KTG': 0.3, '45KTCRF01_AE': 0.0, '45KTCRF02_AG': 0.0, '45KTOther': 0.1, '45KTAll': 0.0, '45KNA': 0.0, '45KNB': 0.1, '45KNC': 0.1, '45KND': 0.0, '45KNF': 0.1, '45KNG': 0.3, '45KNCRF01_AE': 0.0, '45KNCRF02_AG': 0.0, '45KNOther': 0.2, '45KNAll': 0.1, '45KAA': 0.0, '45KAB': 0.0, '45KAC': 0.0, '45KAD': 0.4, '45KAF': 0.1, '45KAG': 0.0, '45KACRF01_AE': 0.0, '45KACRF02_AG': 0.0, '45KAOther': 0.0, '45KAAll': 0.0, '45KEA': 0.0, '45KEB': 0.0, '45KEC': 0.0, '45KED': 0.0, '45KEF': 0.0, '45KEG': 0.0, '45KECRF01_AE': 0.1, '45KECRF02_AG': 0.0, '45KEOther': 0.0, '45KEAll': 0.0, '45KMA': 0.0, '45KMB': 0.0, '45KMC': 0.0, '45KMD': 0.0, '45KMF': 0.0, '45KMG': 0.0, '45KMCRF01_AE': 0.0, '45KMCRF02_AG': 0.0, '45KMOther': 0.0, '45KMAll': 0.0, '45KLA': 0.0, '45KLB': 0.0, '45KLC': 0.0, '45KLD': 0.0, '45KLF': 0.0, '45KLG': 0.0, '45KLCRF01_AE': 0.0, '45KLCRF02_AG': 0.0, '45KLOther': 0.0, '45KLAll': 0.0, '45KPA': 0.0, '45KPB': 0.0, '45KPC': 0.0, '45KPD': 0.0, '45KPF': 0.0, '45KPG': 0.0, '45KPCRF01_AE': 0.0, '45KPCRF02_AG': 0.0, '45KPOther': 0.0, '45KPAll': 0.0, '45KGA': 0.0, '45KGB': 0.0, '45KGC': 0.0, '45KGD': 0.0, '45KGF': 0.0, '45KGG': 0.0, '45KGCRF01_AE': 0.0, '45KGCRF02_AG': 0.0, '45KGOther': 0.0, '45KGAll': 0.0, '46MIA': 11.6, '46MIB': 23.2, '46MIC': 9.5, '46MID': 14.6, '46MIF': 15.6, '46MIG': 17.5, '46MICRF01_AE': 8.0, '46MICRF02_AG': 18.1, '46MIOther': 11.2, '46MIAll': 19.4, '46MLA': 2.1, '46MLB': 9.1, '46MLC': 1.7, '46MLD': 5.9, '46MLF': 11.8, '46MLG': 5.8, '46MLCRF01_AE': 1.5, '46MLCRF02_AG': 3.8, '46MLOther': 8.2, '46MLAll': 7.9, '46MVA': 0.0, '46MVB': 0.3, '46MVC': 0.3, '46MVD': 0.0, '46MVF': 0.3, '46MVG': 0.8, '46MVCRF01_AE': 0.1, '46MVCRF02_AG': 0.4, '46MVOther': 0.2, '46MVAll': 0.3, '46MTA': 0.0, '46MTB': 0.0, '46MTC': 0.1, '46MTD': 0.0, '46MTF': 0.0, '46MTG': 0.0, '46MTCRF01_AE': 0.0, '46MTCRF02_AG': 0.0, '46MTOther': 0.0, '46MTAll': 0.0, '46MKA': 0.0, '46MKB': 0.0, '46MKC': 0.0, '46MKD': 0.0, '46MKF': 0.0, '46MKG': 0.0, '46MKCRF01_AE': 0.0, '46MKCRF02_AG': 0.0, '46MKOther': 0.0, '46MKAll': 0.0, '46MRA': 0.0, '46MRB': 0.0, '46MRC': 0.0, '46MRD': 0.0, '46MRF': 0.0, '46MRG': 0.0, '46MRCRF01_AE': 0.0, '46MRCRF02_AG': 0.0, '46MROther': 0.0, '46MRAll': 0.0, '47IVA': 1.5, '47IVB': 5.2, '47IVC': 0.3, '47IVD': 2.5, '47IVF': 3.8, '47IVG': 2.9, '47IVCRF01_AE': 1.2, '47IVCRF02_AG': 3.4, '47IVOther': 2.3, '47IVAll': 4.1, '47IAA': 0.3, '47IAB': 0.5, '47IAC': 0.5, '47IAD': 0.8, '47IAF': 0.7, '47IAG': 0.5, '47IACRF01_AE': 0.4, '47IACRF02_AG': 0.0, '47IAOther': 0.4, '47IAAll': 0.5, '47IMA': 0.0, '47IMB': 0.0, '47IMC': 0.1, '47IMD': 0.0, '47IMF': 0.1, '47IMG': 0.0, '47IMCRF01_AE': 0.1, '47IMCRF02_AG': 0.0, '47IMOther': 0.0, '47IMAll': 0.0, '47IKA': 0.0, '47IKB': 0.0, '47IKC': 0.0, '47IKD': 0.0, '47IKF': 0.0, '47IKG': 0.0, '47IKCRF01_AE': 0.0, '47IKCRF02_AG': 0.0, '47IKOther': 0.0, '47IKAll': 0.0, '47ILA': 0.0, '47ILB': 0.0, '47ILC': 0.0, '47ILD': 0.0, '47ILF': 0.0, '47ILG': 0.0, '47ILCRF01_AE': 0.0, '47ILCRF02_AG': 0.0, '47ILOther': 0.0, '47ILAll': 0.0, '47ITA': 0.0, '47ITB': 0.0, '47ITC': 0.0, '47ITD': 0.0, '47ITF': 0.0, '47ITG': 0.0, '47ITCRF01_AE': 0.0, '47ITCRF02_AG': 0.0, '47ITOther': 0.0, '47ITAll': 0.0, '47IRA': 0.0, '47IRB': 0.0, '47IRC': 0.0, '47IRD': 0.0, '47IRF': 0.0, '47IRG': 0.0, '47IRCRF01_AE': 0.0, '47IRCRF02_AG': 0.0, '47IROther': 0.0, '47IRAll': 0.0, '48GVA': 1.2, '48GVB': 3.3, '48GVC': 0.7, '48GVD': 0.4, '48GVF': 2.2, '48GVG': 2.3, '48GVCRF01_AE': 4.2, '48GVCRF02_AG': 1.1, '48GVOther': 2.2, '48GVAll': 2.8, '48GAA': 0.0, '48GAB': 0.4, '48GAC': 0.2, '48GAD': 0.4, '48GAF': 0.8, '48GAG': 0.3, '48GACRF01_AE': 0.0, '48GACRF02_AG': 0.8, '48GAOther': 0.3, '48GAAll': 0.4, '48GMA': 0.3, '48GMB': 0.5, '48GMC': 0.0, '48GMD': 0.0, '48GMF': 0.7, '48GMG': 0.5, '48GMCRF01_AE': 0.3, '48GMCRF02_AG': 0.0, '48GMOther': 0.4, '48GMAll': 0.4, '48GQA': 0.3, '48GQB': 0.1, '48GQC': 0.0, '48GQD': 0.0, '48GQF': 0.2, '48GQG': 0.5, '48GQCRF01_AE': 0.0, '48GQCRF02_AG': 0.0, '48GQOther': 0.2, '48GQAll': 0.1, '48GRA': 0.6, '48GRB': 0.0, '48GRC': 0.0, '48GRD': 0.0, '48GRF': 0.1, '48GRG': 0.0, '48GRCRF01_AE': 0.0, '48GRCRF02_AG': 0.4, '48GROther': 0.1, '48GRAll': 0.0, '48GTA': 0.3, '48GTB': 0.0, '48GTC': 0.0, '48GTD': 0.0, '48GTF': 0.1, '48GTG': 0.0, '48GTCRF01_AE': 0.0, '48GTCRF02_AG': 0.4, '48GTOther': 0.0, '48GTAll': 0.0, '48GEA': 0.0, '48GEB': 0.1, '48GEC': 0.1, '48GED': 0.0, '48GEF': 0.1, '48GEG': 0.0, '48GECRF01_AE': 0.0, '48GECRF02_AG': 0.0, '48GEOther': 0.2, '48GEAll': 0.1, '48GSA': 0.0, '48GSB': 0.1, '48GSC': 0.0, '48GSD': 0.0, '48GSF': 0.3, '48GSG': 0.0, '48GSCRF01_AE': 0.0, '48GSCRF02_AG': 0.0, '48GSOther': 0.0, '48GSAll': 0.1, '48GLA': 0.0, '48GLB': 0.1, '48GLC': 0.0, '48GLD': 0.0, '48GLF': 0.0, '48GLG': 0.0, '48GLCRF01_AE': 0.0, '48GLCRF02_AG': 0.0, '48GLOther': 0.0, '48GLAll': 0.1, '48GWA': 0.0, '48GWB': 0.0, '48GWC': 0.0, '48GWD': 0.0, '48GWF': 0.1, '48GWG': 0.0, '48GWCRF01_AE': 0.0, '48GWCRF02_AG': 0.0, '48GWOther': 0.1, '48GWAll': 0.0, '48GIA': 0.0, '48GIB': 0.0, '48GIC': 0.0, '48GID': 0.0, '48GIF': 0.0, '48GIG': 0.0, '48GICRF01_AE': 0.1, '48GICRF02_AG': 0.0, '48GIOther': 0.0, '48GIAll': 0.0, '48GCA': 0.0, '48GCB': 0.0, '48GCC': 0.0, '48GCD': 0.0, '48GCF': 0.0, '48GCG': 0.0, '48GCCRF01_AE': 0.0, '48GCCRF02_AG': 0.0, '48GCOther': 0.0, '48GCAll': 0.0, '49GEA': 0.0, '49GEB': 0.1, '49GEC': 0.0, '49GED': 0.0, '49GEF': 0.1, '49GEG': 0.0, '49GECRF01_AE': 0.0, '49GECRF02_AG': 0.0, '49GEOther': 0.1, '49GEAll': 0.0, '49GRA': 0.0, '49GRB': 0.0, '49GRC': 0.0, '49GRD': 0.0, '49GRF': 0.0, '49GRG': 0.0, '49GRCRF01_AE': 0.0, '49GRCRF02_AG': 0.0, '49GROther': 0.1, '49GRAll': 0.0, '49GKA': 0.0, '49GKB': 0.0, '49GKC': 0.0, '49GKD': 0.0, '49GKF': 0.0, '49GKG': 0.0, '49GKCRF01_AE': 0.0, '49GKCRF02_AG': 0.0, '49GKOther': 0.0, '49GKAll': 0.0, '49GVA': 0.0, '49GVB': 0.0, '49GVC': 0.0, '49GVD': 0.0, '49GVF': 0.0, '49GVG': 0.0, '49GVCRF01_AE': 0.0, '49GVCRF02_AG': 0.0, '49GVOther': 0.0, '49GVAll': 0.0, '49GAA': 0.0, '49GAB': 0.0, '49GAC': 0.0, '49GAD': 0.0, '49GAF': 0.0, '49GAG': 0.0, '49GACRF01_AE': 0.0, '49GACRF02_AG': 0.0, '49GAOther': 0.0, '49GAAll': 0.0, '50ILA': 0.3, '50ILB': 2.2, '50ILC': 2.1, '50ILD': 1.7, '50ILF': 5.2, '50ILG': 0.0, '50ILCRF01_AE': 0.3, '50ILCRF02_AG': 0.0, '50ILOther': 1.1, '50ILAll': 2.1, '50IVA': 0.0, '50IVB': 2.0, '50IVC': 0.6, '50IVD': 2.1, '50IVF': 3.0, '50IVG': 2.3, '50IVCRF01_AE': 0.3, '50IVCRF02_AG': 0.0, '50IVOther': 1.3, '50IVAll': 1.7, '50IMA': 0.0, '50IMB': 0.0, '50IMC': 0.0, '50IMD': 0.4, '50IMF': 0.1, '50IMG': 0.0, '50IMCRF01_AE': 0.0, '50IMCRF02_AG': 0.0, '50IMOther': 0.0, '50IMAll': 0.0, '50INA': 0.0, '50INB': 0.0, '50INC': 0.0, '50IND': 0.0, '50INF': 0.0, '50ING': 0.0, '50INCRF01_AE': 0.1, '50INCRF02_AG': 0.0, '50INOther': 0.1, '50INAll': 0.0, '50IFA': 0.0, '50IFB': 0.0, '50IFC': 0.0, '50IFD': 0.0, '50IFF': 0.0, '50IFG': 0.0, '50IFCRF01_AE': 0.0, '50IFCRF02_AG': 0.0, '50IFOther': 0.0, '50IFAll': 0.0, '50ISA': 0.0, '50ISB': 0.0, '50ISC': 0.0, '50ISD': 0.0, '50ISF': 0.0, '50ISG': 0.0, '50ISCRF01_AE': 0.0, '50ISCRF02_AG': 0.0, '50ISOther': 0.1, '50ISAll': 0.0, '50ITA': 0.0, '50ITB': 0.0, '50ITC': 0.0, '50ITD': 0.0, '50ITF': 0.0, '50ITG': 0.0, '50ITCRF01_AE': 0.0, '50ITCRF02_AG': 0.0, '50ITOther': 0.1, '50ITAll': 0.0, '50IKA': 0.0, '50IKB': 0.0, '50IKC': 0.0, '50IKD': 0.0, '50IKF': 0.0, '50IKG': 0.0, '50IKCRF01_AE': 0.0, '50IKCRF02_AG': 0.0, '50IKOther': 0.0, '50IKAll': 0.0, '51GAA': 0.0, '51GAB': 0.4, '51GAC': 0.1, '51GAD': 0.0, '51GAF': 0.1, '51GAG': 0.0, '51GACRF01_AE': 0.0, '51GACRF02_AG': 0.0, '51GAOther': 0.1, '51GAAll': 0.3, '51GEA': 0.0, '51GEB': 0.0, '51GEC': 0.0, '51GED': 0.0, '51GEF': 0.0, '51GEG': 0.0, '51GECRF01_AE': 0.0, '51GECRF02_AG': 0.0, '51GEOther': 0.1, '51GEAll': 0.0, '51GKA': 0.0, '51GKB': 0.0, '51GKC': 0.0, '51GKD': 0.0, '51GKF': 0.0, '51GKG': 0.0, '51GKCRF01_AE': 0.0, '51GKCRF02_AG': 0.0, '51GKOther': 0.0, '51GKAll': 0.0, '51GVA': 0.0, '51GVB': 0.0, '51GVC': 0.0, '51GVD': 0.0, '51GVF': 0.0, '51GVG': 0.0, '51GVCRF01_AE': 0.0, '51GVCRF02_AG': 0.0, '51GVOther': 0.0, '51GVAll': 0.0, '51GWA': 0.0, '51GWB': 0.0, '51GWC': 0.0, '51GWD': 0.0, '51GWF': 0.0, '51GWG': 0.0, '51GWCRF01_AE': 0.0, '51GWCRF02_AG': 0.0, '51GWOther': 0.0, '51GWAll': 0.0, '51GRA': 0.0, '51GRB': 0.0, '51GRC': 0.0, '51GRD': 0.0, '51GRF': 0.0, '51GRG': 0.0, '51GRCRF01_AE': 0.0, '51GRCRF02_AG': 0.0, '51GROther': 0.0, '51GRAll': 0.0, '52GSA': 0.0, '52GSB': 0.0, '52GSC': 0.0, '52GSD': 0.4, '52GSF': 0.0, '52GSG': 0.0, '52GSCRF01_AE': 0.0, '52GSCRF02_AG': 0.0, '52GSOther': 0.0, '52GSAll': 0.0, '52GEA': 0.0, '52GEB': 0.0, '52GEC': 0.0, '52GED': 0.0, '52GEF': 0.0, '52GEG': 0.0, '52GECRF01_AE': 0.0, '52GECRF02_AG': 0.0, '52GEOther': 0.0, '52GEAll': 0.0, '52GVA': 0.0, '52GVB': 0.0, '52GVC': 0.0, '52GVD': 0.0, '52GVF': 0.0, '52GVG': 0.0, '52GVCRF01_AE': 0.0, '52GVCRF02_AG': 0.0, '52GVOther': 0.0, '52GVAll': 0.0, '52G~A': 0.0, '52G~B': 0.0, '52G~C': 0.0, '52G~D': 0.0, '52G~F': 0.0, '52G~G': 0.0, '52G~CRF01_AE': 0.0, '52G~CRF02_AG': 0.0, '52G~Other': 0.0, '52G~All': 0.0, '52GCA': 0.0, '52GCB': 0.0, '52GCC': 0.0, '52GCD': 0.0, '52GCF': 0.0, '52GCG': 0.0, '52GCCRF01_AE': 0.0, '52GCCRF02_AG': 0.0, '52GCOther': 0.0, '52GCAll': 0.0, '52GAA': 0.0, '52GAB': 0.0, '52GAC': 0.0, '52GAD': 0.0, '52GAF': 0.0, '52GAG': 0.0, '52GACRF01_AE': 0.0, '52GACRF02_AG': 0.0, '52GAOther': 0.0, '52GAAll': 0.0, '52GDA': 0.0, '52GDB': 0.0, '52GDC': 0.0, '52GDD': 0.0, '52GDF': 0.0, '52GDG': 0.0, '52GDCRF01_AE': 0.0, '52GDCRF02_AG': 0.0, '52GDOther': 0.0, '52GDAll': 0.0, '52GRA': 0.0, '52GRB': 0.0, '52GRC': 0.0, '52GRD': 0.0, '52GRF': 0.0, '52GRG': 0.0, '52GRCRF01_AE': 0.0, '52GRCRF02_AG': 0.0, '52GROther': 0.0, '52GRAll': 0.0, '53FLA': 3.9, '53FLB': 5.8, '53FLC': 0.8, '53FLD': 5.9, '53FLF': 10.9, '53FLG': 8.9, '53FLCRF01_AE': 2.8, '53FLCRF02_AG': 2.7, '53FLOther': 6.9, '53FLAll': 5.5, '53FYA': 0.3, '53FYB': 0.4, '53FYC': 0.1, '53FYD': 0.0, '53FYF': 0.3, '53FYG': 0.0, '53FYCRF01_AE': 0.1, '53FYCRF02_AG': 0.4, '53FYOther': 0.1, '53FYAll': 0.3, '53FIA': 0.3, '53FIB': 0.1, '53FIC': 0.0, '53FID': 0.0, '53FIF': 0.5, '53FIG': 0.0, '53FICRF01_AE': 0.1, '53FICRF02_AG': 0.0, '53FIOther': 0.1, '53FIAll': 0.1, '53FCA': 0.0, '53FCB': 0.0, '53FCC': 0.0, '53FCD': 0.0, '53FCF': 0.0, '53FCG': 0.0, '53FCCRF01_AE': 0.0, '53FCCRF02_AG': 0.0, '53FCOther': 0.1, '53FCAll': 0.0, '53FSA': 0.0, '53FSB': 0.0, '53FSC': 0.0, '53FSD': 0.0, '53FSF': 0.0, '53FSG': 0.0, '53FSCRF01_AE': 0.0, '53FSCRF02_AG': 0.0, '53FSOther': 0.0, '53FSAll': 0.0, '53FWA': 0.0, '53FWB': 0.1, '53FWC': 0.0, '53FWD': 0.0, '53FWF': 0.0, '53FWG': 0.0, '53FWCRF01_AE': 0.0, '53FWCRF02_AG': 0.0, '53FWOther': 0.0, '53FWAll': 0.0, '53FVA': 0.0, '53FVB': 0.0, '53FVC': 0.0, '53FVD': 0.0, '53FVF': 0.0, '53FVG': 0.0, '53FVCRF01_AE': 0.0, '53FVCRF02_AG': 0.0, '53FVOther': 0.0, '53FVAll': 0.0, '53F~A': 0.0, '53F~B': 0.0, '53F~C': 0.0, '53F~D': 0.0, '53F~F': 0.0, '53F~G': 0.0, '53F~CRF01_AE': 0.0, '53F~CRF02_AG': 0.0, '53F~Other': 0.0, '53F~All': 0.0, '54IVA': 9.5, '54IVB': 22.8, '54IVC': 11.4, '54IVD': 20.1, '54IVF': 29.0, '54IVG': 33.5, '54IVCRF01_AE': 8.8, '54IVCRF02_AG': 13.8, '54IVOther': 25.7, '54IVAll': 21.5, '54ILA': 0.9, '54ILB': 3.3, '54ILC': 0.7, '54ILD': 3.3, '54ILF': 1.5, '54ILG': 0.5, '54ILCRF01_AE': 0.0, '54ILCRF02_AG': 0.8, '54ILOther': 1.2, '54ILAll': 2.6, '54IMA': 0.3, '54IMB': 2.6, '54IMC': 0.2, '54IMD': 0.8, '54IMF': 0.7, '54IMG': 0.3, '54IMCRF01_AE': 0.1, '54IMCRF02_AG': 0.0, '54IMOther': 0.6, '54IMAll': 1.8, '54IAA': 0.0, '54IAB': 1.1, '54IAC': 0.2, '54IAD': 0.4, '54IAF': 1.2, '54IAG': 0.8, '54IACRF01_AE': 1.0, '54IACRF02_AG': 0.0, '54IAOther': 1.0, '54IAAll': 0.9, '54ITA': 1.2, '54ITB': 0.8, '54ITC': 0.2, '54ITD': 0.0, '54ITF': 0.4, '54ITG': 0.3, '54ITCRF01_AE': 0.1, '54ITCRF02_AG': 0.0, '54ITOther': 0.2, '54ITAll': 0.6, '54ISA': 0.3, '54ISB': 0.7, '54ISC': 0.0, '54ISD': 0.4, '54ISF': 0.4, '54ISG': 0.0, '54ISCRF01_AE': 0.1, '54ISCRF02_AG': 0.0, '54ISOther': 0.3, '54ISAll': 0.5, '54ICA': 0.0, '54ICB': 0.0, '54ICC': 0.1, '54ICD': 0.0, '54ICF': 0.0, '54ICG': 0.0, '54ICCRF01_AE': 0.1, '54ICCRF02_AG': 0.0, '54ICOther': 0.1, '54ICAll': 0.0, '54IFA': 0.3, '54IFB': 0.0, '54IFC': 0.0, '54IFD': 0.0, '54IFF': 0.0, '54IFG': 0.0, '54IFCRF01_AE': 0.0, '54IFCRF02_AG': 0.0, '54IFOther': 0.0, '54IFAll': 0.0, '54INA': 0.0, '54INB': 0.0, '54INC': 0.0, '54IND': 0.0, '54INF': 0.0, '54ING': 0.0, '54INCRF01_AE': 0.0, '54INCRF02_AG': 0.0, '54INOther': 0.0, '54INAll': 0.0, '54IKA': 0.0, '54IKB': 0.0, '54IKC': 0.0, '54IKD': 0.0, '54IKF': 0.0, '54IKG': 0.0, '54IKCRF01_AE': 0.0, '54IKCRF02_AG': 0.0, '54IKOther': 0.0, '54IKAll': 0.0, '54IPA': 0.0, '54IPB': 0.0, '54IPC': 0.0, '54IPD': 0.0, '54IPF': 0.0, '54IPG': 0.0, '54IPCRF01_AE': 0.0, '54IPCRF02_AG': 0.0, '54IPOther': 0.0, '54IPAll': 0.0, '55KRA': 2.1, '55KRB': 8.0, '55KRC': 2.1, '55KRD': 4.2, '55KRF': 13.8, '55KRG': 2.3, '55KRCRF01_AE': 3.6, '55KRCRF02_AG': 1.9, '55KROther': 8.3, '55KRAll': 7.3, '55KNA': 0.6, '55KNB': 0.4, '55KNC': 0.3, '55KND': 0.0, '55KNF': 0.6, '55KNG': 0.8, '55KNCRF01_AE': 0.0, '55KNCRF02_AG': 0.0, '55KNOther': 0.3, '55KNAll': 0.3, '55KHA': 0.0, '55KHB': 0.0, '55KHC': 0.1, '55KHD': 0.4, '55KHF': 0.0, '55KHG': 0.3, '55KHCRF01_AE': 0.0, '55KHCRF02_AG': 0.0, '55KHOther': 0.1, '55KHAll': 0.1, '55KTA': 0.0, '55KTB': 0.0, '55KTC': 0.0, '55KTD': 0.4, '55KTF': 0.1, '55KTG': 0.0, '55KTCRF01_AE': 0.0, '55KTCRF02_AG': 0.0, '55KTOther': 0.2, '55KTAll': 0.0, '55KIA': 0.0, '55KIB': 0.0, '55KIC': 0.0, '55KID': 0.4, '55KIF': 0.0, '55KIG': 0.0, '55KICRF01_AE': 0.0, '55KICRF02_AG': 0.0, '55KIOther': 0.0, '55KIAll': 0.0, '55KQA': 0.0, '55KQB': 0.1, '55KQC': 0.0, '55KQD': 0.0, '55KQF': 0.1, '55KQG': 0.0, '55KQCRF01_AE': 0.0, '55KQCRF02_AG': 0.0, '55KQOther': 0.0, '55KQAll': 0.0, '55KFA': 0.0, '55KFB': 0.0, '55KFC': 0.0, '55KFD': 0.0, '55KFF': 0.0, '55KFG': 0.0, '55KFCRF01_AE': 0.0, '55KFCRF02_AG': 0.0, '55KFOther': 0.0, '55KFAll': 0.0, '55KEA': 0.0, '55KEB': 0.0, '55KEC': 0.0, '55KED': 0.0, '55KEF': 0.1, '55KEG': 0.0, '55KECRF01_AE': 0.0, '55KECRF02_AG': 0.0, '55KEOther': 0.0, '55KEAll': 0.0, '55KMA': 0.0, '55KMB': 0.0, '55KMC': 0.0, '55KMD': 0.0, '55KMF': 0.0, '55KMG': 0.0, '55KMCRF01_AE': 0.0, '55KMCRF02_AG': 0.0, '55KMOther': 0.0, '55KMAll': 0.0, '55KLA': 0.0, '55KLB': 0.0, '55KLC': 0.0, '55KLD': 0.0, '55KLF': 0.0, '55KLG': 0.0, '55KLCRF01_AE': 0.0, '55KLCRF02_AG': 0.0, '55KLOther': 0.0, '55KLAll': 0.0, '56VAA': 0.0, '56VAB': 0.0, '56VAC': 0.0, '56VAD': 0.4, '56VAF': 0.3, '56VAG': 0.3, '56VACRF01_AE': 0.0, '56VACRF02_AG': 0.0, '56VAOther': 0.1, '56VAAll': 0.1, '56VIA': 0.0, '56VIB': 0.0, '56VIC': 0.0, '56VID': 0.0, '56VIF': 0.1, '56VIG': 0.0, '56VICRF01_AE': 0.0, '56VICRF02_AG': 0.0, '56VIOther': 0.1, '56VIAll': 0.0, '56VLA': 0.0, '56VLB': 0.1, '56VLC': 0.0, '56VLD': 0.0, '56VLF': 0.0, '56VLG': 0.0, '56VLCRF01_AE': 0.0, '56VLCRF02_AG': 0.0, '56VLOther': 0.0, '56VLAll': 0.0, '56VGA': 0.0, '56VGB': 0.0, '56VGC': 0.0, '56VGD': 0.0, '56VGF': 0.0, '56VGG': 0.0, '56VGCRF01_AE': 0.0, '56VGCRF02_AG': 0.0, '56VGOther': 0.0, '56VGAll': 0.0, '56VTA': 0.0, '56VTB': 0.0, '56VTC': 0.0, '56VTD': 0.0, '56VTF': 0.0, '56VTG': 0.0, '56VTCRF01_AE': 0.0, '56VTCRF02_AG': 0.0, '56VTOther': 0.0, '56VTAll': 0.0, '56VKA': 0.0, '56VKB': 0.0, '56VKC': 0.0, '56VKD': 0.0, '56VKF': 0.0, '56VKG': 0.0, '56VKCRF01_AE': 0.0, '56VKCRF02_AG': 0.0, '56VKOther': 0.0, '56VKAll': 0.0, '56VEA': 0.0, '56VEB': 0.0, '56VEC': 0.0, '56VED': 0.0, '56VEF': 0.0, '56VEG': 0.0, '56VECRF01_AE': 0.0, '56VECRF02_AG': 0.0, '56VEOther': 0.0, '56VEAll': 0.0, '56VMA': 0.0, '56VMB': 0.0, '56VMC': 0.0, '56VMD': 0.0, '56VMF': 0.0, '56VMG': 0.0, '56VMCRF01_AE': 0.0, '56VMCRF02_AG': 0.0, '56VMOther': 0.0, '56VMAll': 0.0, '56VRA': 0.0, '56VRB': 0.0, '56VRC': 0.0, '56VRD': 0.0, '56VRF': 0.0, '56VRG': 0.0, '56VRCRF01_AE': 0.0, '56VRCRF02_AG': 0.0, '56VROther': 0.0, '56VRAll': 0.0, '57RKA': 40.9, '57RKB': 16.2, '57RKC': 8.6, '57RKD': 13.4, '57RKF': 93.2, '57RKG': 13.0, '57RKCRF01_AE': 10.5, '57RKCRF02_AG': 3.4, '57RKOther': 85.9, '57RKAll': 25.3, '57RGA': 0.3, '57RGB': 0.0, '57RGC': 0.0, '57RGD': 0.0, '57RGF': 0.0, '57RGG': 0.3, '57RGCRF01_AE': 0.0, '57RGCRF02_AG': 0.0, '57RGOther': 0.0, '57RGAll': 0.0, '57RSA': 0.0, '57RSB': 0.0, '57RSC': 0.0, '57RSD': 0.4, '57RSF': 0.0, '57RSG': 0.0, '57RSCRF01_AE': 0.0, '57RSCRF02_AG': 0.0, '57RSOther': 0.1, '57RSAll': 0.0, '57RQA': 0.3, '57RQB': 0.0, '57RQC': 0.0, '57RQD': 0.0, '57RQF': 0.0, '57RQG': 0.0, '57RQCRF01_AE': 0.0, '57RQCRF02_AG': 0.0, '57RQOther': 0.0, '57RQAll': 0.0, '57RNA': 0.3, '57RNB': 0.0, '57RNC': 0.0, '57RND': 0.0, '57RNF': 0.0, '57RNG': 0.0, '57RNCRF01_AE': 0.0, '57RNCRF02_AG': 0.0, '57RNOther': 0.0, '57RNAll': 0.0, '57REA': 0.0, '57REB': 0.0, '57REC': 0.0, '57RED': 0.0, '57REF': 0.1, '57REG': 0.0, '57RECRF01_AE': 0.0, '57RECRF02_AG': 0.0, '57REOther': 0.0, '57REAll': 0.0, '57RTA': 0.0, '57RTB': 0.0, '57RTC': 0.0, '57RTD': 0.0, '57RTF': 0.0, '57RTG': 0.0, '57RTCRF01_AE': 0.0, '57RTCRF02_AG': 0.0, '57RTOther': 0.0, '57RTAll': 0.0, '57RPA': 0.0, '57RPB': 0.0, '57RPC': 0.0, '57RPD': 0.0, '57RPF': 0.0, '57RPG': 0.0, '57RPCRF01_AE': 0.0, '57RPCRF02_AG': 0.0, '57RPOther': 0.0, '57RPAll': 0.0, '57RIA': 0.0, '57RIB': 0.0, '57RIC': 0.0, '57RID': 0.0, '57RIF': 0.0, '57RIG': 0.0, '57RICRF01_AE': 0.0, '57RICRF02_AG': 0.0, '57RIOther': 0.0, '57RIAll': 0.0, '58QEA': 2.4, '58QEB': 7.6, '58QEC': 2.5, '58QED': 6.3, '58QEF': 5.4, '58QEG': 2.9, '58QECRF01_AE': 2.5, '58QECRF02_AG': 2.7, '58QEOther': 4.2, '58QEAll': 6.3, '58QRA': 0.3, '58QRB': 0.0, '58QRC': 0.0, '58QRD': 0.0, '58QRF': 0.0, '58QRG': 0.3, '58QRCRF01_AE': 0.0, '58QRCRF02_AG': 0.0, '58QROther': 0.0, '58QRAll': 0.0, '58QKA': 0.0, '58QKB': 0.0, '58QKC': 0.0, '58QKD': 0.0, '58QKF': 0.1, '58QKG': 0.0, '58QKCRF01_AE': 0.0, '58QKCRF02_AG': 0.0, '58QKOther': 0.0, '58QKAll': 0.0, '58QPA': 0.3, '58QPB': 0.0, '58QPC': 0.0, '58QPD': 0.0, '58QPF': 0.0, '58QPG': 0.0, '58QPCRF01_AE': 0.0, '58QPCRF02_AG': 0.0, '58QPOther': 0.0, '58QPAll': 0.0, '58QHA': 0.0, '58QHB': 0.0, '58QHC': 0.0, '58QHD': 0.0, '58QHF': 0.0, '58QHG': 0.0, '58QHCRF01_AE': 0.0, '58QHCRF02_AG': 0.0, '58QHOther': 0.0, '58QHAll': 0.0, '58QTA': 0.0, '58QTB': 0.0, '58QTC': 0.0, '58QTD': 0.0, '58QTF': 0.0, '58QTG': 0.0, '58QTCRF01_AE': 0.0, '58QTCRF02_AG': 0.0, '58QTOther': 0.0, '58QTAll': 0.0, '58QLA': 0.0, '58QLB': 0.0, '58QLC': 0.0, '58QLD': 0.0, '58QLF': 0.0, '58QLG': 0.0, '58QLCRF01_AE': 0.0, '58QLCRF02_AG': 0.0, '58QLOther': 0.0, '58QLAll': 0.0, '59YFA': 0.3, '59YFB': 0.1, '59YFC': 0.1, '59YFD': 0.0, '59YFF': 0.1, '59YFG': 0.0, '59YFCRF01_AE': 0.0, '59YFCRF02_AG': 0.0, '59YFOther': 0.1, '59YFAll': 0.1, '59YHA': 0.0, '59YHB': 0.1, '59YHC': 0.0, '59YHD': 0.0, '59YHF': 0.2, '59YHG': 0.0, '59YHCRF01_AE': 0.0, '59YHCRF02_AG': 0.0, '59YHOther': 0.2, '59YHAll': 0.1, '59YCA': 0.0, '59YCB': 0.0, '59YCC': 0.0, '59YCD': 0.0, '59YCF': 0.0, '59YCG': 0.0, '59YCCRF01_AE': 0.0, '59YCCRF02_AG': 0.4, '59YCOther': 0.0, '59YCAll': 0.0, '59YSA': 0.0, '59YSB': 0.0, '59YSC': 0.0, '59YSD': 0.0, '59YSF': 0.0, '59YSG': 0.0, '59YSCRF01_AE': 0.0, '59YSCRF02_AG': 0.0, '59YSOther': 0.1, '59YSAll': 0.0, '59YNA': 0.0, '59YNB': 0.0, '59YNC': 0.0, '59YND': 0.0, '59YNF': 0.0, '59YNG': 0.0, '59YNCRF01_AE': 0.0, '59YNCRF02_AG': 0.0, '59YNOther': 0.0, '59YNAll': 0.0, '59YEA': 0.0, '59YEB': 0.0, '59YEC': 0.0, '59YED': 0.0, '59YEF': 0.0, '59YEG': 0.0, '59YECRF01_AE': 0.0, '59YECRF02_AG': 0.0, '59YEOther': 0.0, '59YEAll': 0.0, '59YVA': 0.0, '59YVB': 0.0, '59YVC': 0.0, '59YVD': 0.0, '59YVF': 0.0, '59YVG': 0.0, '59YVCRF01_AE': 0.0, '59YVCRF02_AG': 0.0, '59YVOther': 0.0, '59YVAll': 0.0, '59YMA': 0.0, '59YMB': 0.0, '59YMC': 0.0, '59YMD': 0.0, '59YMF': 0.0, '59YMG': 0.0, '59YMCRF01_AE': 0.0, '59YMCRF02_AG': 0.0, '59YMOther': 0.0, '59YMAll': 0.0, '59YDA': 0.0, '59YDB': 0.0, '59YDC': 0.0, '59YDD': 0.0, '59YDF': 0.0, '59YDG': 0.0, '59YDCRF01_AE': 0.0, '59YDCRF02_AG': 0.0, '59YDOther': 0.0, '59YDAll': 0.0, '59YIA': 0.0, '59YIB': 0.0, '59YIC': 0.0, '59YID': 0.0, '59YIF': 0.0, '59YIG': 0.0, '59YICRF01_AE': 0.0, '59YICRF02_AG': 0.0, '59YIOther': 0.0, '59YIAll': 0.0, '60DEA': 7.1, '60DEB': 11.8, '60DEC': 14.8, '60DED': 14.6, '60DEF': 20.4, '60DEG': 5.7, '60DECRF01_AE': 3.4, '60DECRF02_AG': 0.8, '60DEOther': 15.5, '60DEAll': 12.5, '60DNA': 0.6, '60DNB': 0.2, '60DNC': 0.6, '60DND': 0.0, '60DNF': 0.4, '60DNG': 0.3, '60DNCRF01_AE': 0.0, '60DNCRF02_AG': 0.0, '60DNOther': 0.6, '60DNAll': 0.3, '60DGA': 0.0, '60DGB': 0.0, '60DGC': 0.2, '60DGD': 0.0, '60DGF': 0.1, '60DGG': 0.0, '60DGCRF01_AE': 0.0, '60DGCRF02_AG': 0.0, '60DGOther': 0.0, '60DGAll': 0.0, '60DSA': 0.0, '60DSB': 0.0, '60DSC': 0.0, '60DSD': 0.0, '60DSF': 0.1, '60DSG': 0.0, '60DSCRF01_AE': 0.0, '60DSCRF02_AG': 0.0, '60DSOther': 0.0, '60DSAll': 0.0, '60DQA': 0.0, '60DQB': 0.0, '60DQC': 0.0, '60DQD': 0.0, '60DQF': 0.0, '60DQG': 0.3, '60DQCRF01_AE': 0.0, '60DQCRF02_AG': 0.0, '60DQOther': 0.0, '60DQAll': 0.0, '60DKA': 0.0, '60DKB': 0.0, '60DKC': 0.1, '60DKD': 0.0, '60DKF': 0.1, '60DKG': 0.0, '60DKCRF01_AE': 0.0, '60DKCRF02_AG': 0.0, '60DKOther': 0.0, '60DKAll': 0.0, '60DTA': 0.0, '60DTB': 0.0, '60DTC': 0.0, '60DTD': 0.0, '60DTF': 0.0, '60DTG': 0.0, '60DTCRF01_AE': 0.0, '60DTCRF02_AG': 0.0, '60DTOther': 0.1, '60DTAll': 0.0, '60DHA': 0.0, '60DHB': 0.0, '60DHC': 0.0, '60DHD': 0.0, '60DHF': 0.0, '60DHG': 0.0, '60DHCRF01_AE': 0.1, '60DHCRF02_AG': 0.0, '60DHOther': 0.0, '60DHAll': 0.0, '60DRA': 0.0, '60DRB': 0.0, '60DRC': 0.0, '60DRD': 0.0, '60DRF': 0.0, '60DRG': 0.0, '60DRCRF01_AE': 0.0, '60DRCRF02_AG': 0.0, '60DROther': 0.1, '60DRAll': 0.0, '60DYA': 0.0, '60DYB': 0.0, '60DYC': 0.0, '60DYD': 0.0, '60DYF': 0.0, '60DYG': 0.0, '60DYCRF01_AE': 0.0, '60DYCRF02_AG': 0.0, '60DYOther': 0.0, '60DYAll': 0.0, '60DVA': 0.0, '60DVB': 0.0, '60DVC': 0.0, '60DVD': 0.0, '60DVF': 0.0, '60DVG': 0.0, '60DVCRF01_AE': 0.0, '60DVCRF02_AG': 0.0, '60DVOther': 0.0, '60DVAll': 0.0, '60DAA': 0.0, '60DAB': 0.0, '60DAC': 0.0, '60DAD': 0.0, '60DAF': 0.0, '60DAG': 0.0, '60DACRF01_AE': 0.0, '60DACRF02_AG': 0.0, '60DAOther': 0.0, '60DAAll': 0.0, '61QNA': 0.6, '61QNB': 1.5, '61QNC': 2.2, '61QND': 7.5, '61QNF': 70.0, '61QNG': 1.3, '61QNCRF01_AE': 0.3, '61QNCRF02_AG': 0.0, '61QNOther': 61.4, '61QNAll': 10.3, '61QDA': 0.0, '61QDB': 0.7, '61QDC': 0.8, '61QDD': 0.8, '61QDF': 21.4, '61QDG': 0.3, '61QDCRF01_AE': 0.4, '61QDCRF02_AG': 0.0, '61QDOther': 23.6, '61QDAll': 3.8, '61QEA': 5.6, '61QEB': 3.1, '61QEC': 8.7, '61QED': 3.3, '61QEF': 1.6, '61QEG': 2.1, '61QECRF01_AE': 1.9, '61QECRF02_AG': 0.0, '61QEOther': 1.4, '61QEAll': 3.5, '61QHA': 1.5, '61QHB': 1.4, '61QHC': 2.7, '61QHD': 0.4, '61QHF': 1.2, '61QHG': 4.2, '61QHCRF01_AE': 0.9, '61QHCRF02_AG': 4.6, '61QHOther': 0.8, '61QHAll': 1.6, '61QSA': 0.0, '61QSB': 0.0, '61QSC': 0.0, '61QSD': 0.0, '61QSF': 1.2, '61QSG': 0.0, '61QSCRF01_AE': 0.0, '61QSCRF02_AG': 0.0, '61QSOther': 1.1, '61QSAll': 0.2, '61QGA': 0.3, '61QGB': 0.0, '61QGC': 0.0, '61QGD': 0.0, '61QGF': 0.1, '61QGG': 0.0, '61QGCRF01_AE': 0.0, '61QGCRF02_AG': 0.0, '61QGOther': 0.2, '61QGAll': 0.1, '61QRA': 0.6, '61QRB': 0.3, '61QRC': 0.1, '61QRD': 0.4, '61QRF': 0.0, '61QRG': 0.0, '61QRCRF01_AE': 0.0, '61QRCRF02_AG': 0.4, '61QROther': 0.1, '61QRAll': 0.2, '61QKA': 0.0, '61QKB': 0.1, '61QKC': 0.2, '61QKD': 0.0, '61QKF': 0.1, '61QKG': 0.0, '61QKCRF01_AE': 0.1, '61QKCRF02_AG': 0.0, '61QKOther': 0.1, '61QKAll': 0.1, '61QYA': 0.0, '61QYB': 0.0, '61QYC': 0.0, '61QYD': 0.0, '61QYF': 0.0, '61QYG': 0.0, '61QYCRF01_AE': 0.0, '61QYCRF02_AG': 0.0, '61QYOther': 0.1, '61QYAll': 0.0, '61QLA': 0.0, '61QLB': 0.0, '61QLC': 0.0, '61QLD': 0.0, '61QLF': 0.0, '61QLG': 0.0, '61QLCRF01_AE': 0.0, '61QLCRF02_AG': 0.0, '61QLOther': 0.0, '61QLAll': 0.0, '61QAA': 0.0, '61QAB': 0.0, '61QAC': 0.0, '61QAD': 0.0, '61QAF': 0.0, '61QAG': 0.0, '61QACRF01_AE': 0.0, '61QACRF02_AG': 0.0, '61QAOther': 0.0, '61QAAll': 0.0, '61QPA': 0.0, '61QPB': 0.0, '61QPC': 0.0, '61QPD': 0.0, '61QPF': 0.0, '61QPG': 0.0, '61QPCRF01_AE': 0.0, '61QPCRF02_AG': 0.0, '61QPOther': 0.0, '61QPAll': 0.0, '62IVA': 10.7, '62IVB': 40.0, '62IVC': 11.3, '62IVD': 33.1, '62IVF': 24.7, '62IVG': 13.8, '62IVCRF01_AE': 11.6, '62IVCRF02_AG': 8.4, '62IVOther': 18.2, '62IVAll': 32.1, '62IMA': 0.0, '62IMB': 0.1, '62IMC': 0.0, '62IMD': 0.0, '62IMF': 0.0, '62IMG': 0.0, '62IMCRF01_AE': 0.1, '62IMCRF02_AG': 0.0, '62IMOther': 0.3, '62IMAll': 0.1, '62ITA': 0.3, '62ITB': 0.0, '62ITC': 0.0, '62ITD': 0.0, '62ITF': 0.1, '62ITG': 0.0, '62ITCRF01_AE': 0.1, '62ITCRF02_AG': 0.0, '62ITOther': 0.2, '62ITAll': 0.1, '62ILA': 0.0, '62ILB': 0.0, '62ILC': 0.1, '62ILD': 0.0, '62ILF': 0.1, '62ILG': 0.0, '62ILCRF01_AE': 0.0, '62ILCRF02_AG': 0.0, '62ILOther': 0.2, '62ILAll': 0.0, '62IKA': 0.0, '62IKB': 0.0, '62IKC': 0.0, '62IKD': 0.0, '62IKF': 0.0, '62IKG': 0.0, '62IKCRF01_AE': 0.0, '62IKCRF02_AG': 0.0, '62IKOther': 0.1, '62IKAll': 0.0, '62IRA': 0.0, '62IRB': 0.0, '62IRC': 0.0, '62IRD': 0.0, '62IRF': 0.1, '62IRG': 0.0, '62IRCRF01_AE': 0.0, '62IRCRF02_AG': 0.0, '62IROther': 0.0, '62IRAll': 0.0, '62IGA': 0.0, '62IGB': 0.0, '62IGC': 0.0, '62IGD': 0.0, '62IGF': 0.0, '62IGG': 0.0, '62IGCRF01_AE': 0.0, '62IGCRF02_AG': 0.0, '62IGOther': 0.0, '62IGAll': 0.0, '62IFA': 0.0, '62IFB': 0.0, '62IFC': 0.0, '62IFD': 0.0, '62IFF': 0.0, '62IFG': 0.0, '62IFCRF01_AE': 0.0, '62IFCRF02_AG': 0.0, '62IFOther': 0.0, '62IFAll': 0.0, '62ISA': 0.0, '62ISB': 0.0, '62ISC': 0.0, '62ISD': 0.0, '62ISF': 0.0, '62ISG': 0.0, '62ISCRF01_AE': 0.0, '62ISCRF02_AG': 0.0, '62ISOther': 0.0, '62ISAll': 0.0, '62IQA': 0.0, '62IQB': 0.0, '62IQC': 0.0, '62IQD': 0.0, '62IQF': 0.0, '62IQG': 0.0, '62IQCRF01_AE': 0.0, '62IQCRF02_AG': 0.0, '62IQOther': 0.0, '62IQAll': 0.0, '63LPA': 19.9, '63LPB': 71.5, '63LPC': 41.2, '63LPD': 70.5, '63LPF': 26.5, '63LPG': 29.2, '63LPCRF01_AE': 20.8, '63LPCRF02_AG': 26.5, '63LPOther': 28.8, '63LPAll': 58.8, '63LTA': 6.0, '63LTB': 3.2, '63LTC': 8.7, '63LTD': 2.1, '63LTF': 13.8, '63LTG': 0.5, '63LTCRF01_AE': 4.3, '63LTCRF02_AG': 1.9, '63LTOther': 6.5, '63LTAll': 4.7, '63LSA': 1.8, '63LSB': 2.9, '63LSC': 2.9, '63LSD': 1.7, '63LSF': 4.9, '63LSG': 2.3, '63LSCRF01_AE': 2.2, '63LSCRF02_AG': 3.8, '63LSOther': 4.6, '63LSAll': 3.1, '63LVA': 6.9, '63LVB': 0.6, '63LVC': 8.3, '63LVD': 1.3, '63LVF': 3.4, '63LVG': 1.3, '63LVCRF01_AE': 0.9, '63LVCRF02_AG': 1.9, '63LVOther': 4.3, '63LVAll': 2.0, '63LAA': 1.8, '63LAB': 3.6, '63LAC': 2.8, '63LAD': 0.8, '63LAF': 3.0, '63LAG': 0.3, '63LACRF01_AE': 0.1, '63LACRF02_AG': 0.0, '63LAOther': 4.2, '63LAAll': 3.4, '63LHA': 0.9, '63LHB': 1.8, '63LHC': 1.6, '63LHD': 1.3, '63LHF': 1.8, '63LHG': 1.6, '63LHCRF01_AE': 0.9, '63LHCRF02_AG': 3.5, '63LHOther': 1.3, '63LHAll': 1.7, '63LCA': 0.6, '63LCB': 1.2, '63LCC': 0.7, '63LCD': 0.8, '63LCF': 4.2, '63LCG': 1.0, '63LCCRF01_AE': 2.7, '63LCCRF02_AG': 0.0, '63LCOther': 4.0, '63LCAll': 1.6, '63LQA': 0.6, '63LQB': 1.9, '63LQC': 0.9, '63LQD': 3.4, '63LQF': 0.8, '63LQG': 0.8, '63LQCRF01_AE': 0.6, '63LQCRF02_AG': 0.0, '63LQOther': 1.0, '63LQAll': 1.5, '63LIA': 1.5, '63LIB': 0.1, '63LIC': 1.1, '63LID': 0.0, '63LIF': 0.9, '63LIG': 0.5, '63LICRF01_AE': 0.4, '63LICRF02_AG': 0.8, '63LIOther': 0.7, '63LIAll': 0.3, '63LNA': 0.6, '63LNB': 0.3, '63LNC': 0.2, '63LND': 0.4, '63LNF': 0.3, '63LNG': 0.0, '63LNCRF01_AE': 0.1, '63LNCRF02_AG': 0.0, '63LNOther': 0.6, '63LNAll': 0.3, '63LRA': 0.6, '63LRB': 0.3, '63LRC': 0.3, '63LRD': 1.3, '63LRF': 0.1, '63LRG': 0.0, '63LRCRF01_AE': 0.0, '63LRCRF02_AG': 1.2, '63LROther': 0.3, '63LRAll': 0.3, '63LGA': 0.0, '63LGB': 0.0, '63LGC': 0.0, '63LGD': 0.0, '63LGF': 0.9, '63LGG': 0.0, '63LGCRF01_AE': 0.0, '63LGCRF02_AG': 0.4, '63LGOther': 1.2, '63LGAll': 0.2, '63LMA': 0.3, '63LMB': 0.1, '63LMC': 0.2, '63LMD': 0.0, '63LMF': 0.0, '63LMG': 0.0, '63LMCRF01_AE': 0.0, '63LMCRF02_AG': 0.4, '63LMOther': 0.1, '63LMAll': 0.1, '63LDA': 0.3, '63LDB': 0.1, '63LDC': 0.0, '63LDD': 0.0, '63LDF': 0.6, '63LDG': 0.0, '63LDCRF01_AE': 0.3, '63LDCRF02_AG': 0.0, '63LDOther': 0.5, '63LDAll': 0.2, '63LFA': 0.0, '63LFB': 0.0, '63LFC': 0.0, '63LFD': 0.4, '63LFF': 0.1, '63LFG': 0.0, '63LFCRF01_AE': 0.1, '63LFCRF02_AG': 0.0, '63LFOther': 0.0, '63LFAll': 0.0, '63LEA': 0.0, '63LEB': 0.3, '63LEC': 0.1, '63LED': 0.0, '63LEF': 0.2, '63LEG': 0.0, '63LECRF01_AE': 0.1, '63LECRF02_AG': 0.0, '63LEOther': 0.2, '63LEAll': 0.2, '63LYA': 0.0, '63LYB': 0.0, '63LYC': 0.0, '63LYD': 0.0, '63LYF': 0.0, '63LYG': 0.0, '63LYCRF01_AE': 0.0, '63LYCRF02_AG': 0.0, '63LYOther': 0.0, '63LYAll': 0.0, '63LKA': 0.0, '63LKB': 0.0, '63LKC': 0.0, '63LKD': 0.0, '63LKF': 0.0, '63LKG': 0.0, '63LKCRF01_AE': 0.0, '63LKCRF02_AG': 0.0, '63LKOther': 0.0, '63LKAll': 0.0, '63LWA': 0.0, '63LWB': 0.0, '63LWC': 0.0, '63LWD': 0.0, '63LWF': 0.0, '63LWG': 0.0, '63LWCRF01_AE': 0.0, '63LWCRF02_AG': 0.0, '63LWOther': 0.0, '63LWAll': 0.0, '63L#A': 0.0, '63L#B': 0.0, '63L#C': 0.0, '63L#D': 0.0, '63L#F': 0.0, '63L#G': 0.0, '63L#CRF01_AE': 0.0, '63L#CRF02_AG': 0.0, '63L#Other': 0.0, '63L#All': 0.0, '64IVA': 0.9, '64IVB': 21.8, '64IVC': 1.8, '64IVD': 29.3, '64IVF': 5.2, '64IVG': 2.6, '64IVCRF01_AE': 3.4, '64IVCRF02_AG': 1.9, '64IVOther': 6.5, '64IVAll': 16.0, '64IMA': 0.6, '64IMB': 1.3, '64IMC': 1.2, '64IMD': 0.8, '64IMF': 1.2, '64IMG': 6.0, '64IMCRF01_AE': 0.4, '64IMCRF02_AG': 18.4, '64IMOther': 1.5, '64IMAll': 1.5, '64ILA': 0.9, '64ILB': 2.7, '64ILC': 1.7, '64ILD': 5.0, '64ILF': 4.1, '64ILG': 1.0, '64ILCRF01_AE': 0.6, '64ILCRF02_AG': 3.8, '64ILOther': 2.7, '64ILAll': 2.6, '64ITA': 0.0, '64ITB': 0.1, '64ITC': 0.0, '64ITD': 0.0, '64ITF': 0.0, '64ITG': 0.0, '64ITCRF01_AE': 0.0, '64ITCRF02_AG': 0.0, '64ITOther': 0.0, '64ITAll': 0.0, '64IRA': 0.0, '64IRB': 0.0, '64IRC': 0.0, '64IRD': 0.0, '64IRF': 0.0, '64IRG': 0.0, '64IRCRF01_AE': 0.0, '64IRCRF02_AG': 0.0, '64IROther': 0.1, '64IRAll': 0.0, '64IFA': 0.0, '64IFB': 0.0, '64IFC': 0.0, '64IFD': 0.0, '64IFF': 0.0, '64IFG': 0.0, '64IFCRF01_AE': 0.0, '64IFCRF02_AG': 0.0, '64IFOther': 0.0, '64IFAll': 0.0, '64IKA': 0.0, '64IKB': 0.0, '64IKC': 0.0, '64IKD': 0.0, '64IKF': 0.0, '64IKG': 0.0, '64IKCRF01_AE': 0.0, '64IKCRF02_AG': 0.0, '64IKOther': 0.0, '64IKAll': 0.0, '64IPA': 0.0, '64IPB': 0.0, '64IPC': 0.0, '64IPD': 0.0, '64IPF': 0.0, '64IPG': 0.0, '64IPCRF01_AE': 0.0, '64IPCRF02_AG': 0.0, '64IPOther': 0.0, '64IPAll': 0.0, '65EDA': 2.4, '65EDB': 3.0, '65EDC': 1.0, '65EDD': 5.9, '65EDF': 21.5, '65EDG': 0.0, '65EDCRF01_AE': 0.6, '65EDCRF02_AG': 0.4, '65EDOther': 10.8, '65EDAll': 4.4, '65EKA': 0.3, '65EKB': 0.1, '65EKC': 0.0, '65EKD': 0.0, '65EKF': 0.0, '65EKG': 0.0, '65EKCRF01_AE': 0.0, '65EKCRF02_AG': 0.4, '65EKOther': 0.4, '65EKAll': 0.1, '65ENA': 0.0, '65ENB': 0.0, '65ENC': 0.0, '65END': 0.0, '65ENF': 0.1, '65ENG': 0.0, '65ENCRF01_AE': 0.0, '65ENCRF02_AG': 0.0, '65ENOther': 0.0, '65ENAll': 0.0, '65EQA': 0.0, '65EQB': 0.0, '65EQC': 0.0, '65EQD': 0.0, '65EQF': 0.0, '65EQG': 0.0, '65EQCRF01_AE': 0.0, '65EQCRF02_AG': 0.4, '65EQOther': 0.0, '65EQAll': 0.0, '65ERA': 0.0, '65ERB': 0.0, '65ERC': 0.0, '65ERD': 0.0, '65ERF': 0.1, '65ERG': 0.0, '65ERCRF01_AE': 0.0, '65ERCRF02_AG': 0.0, '65EROther': 0.0, '65ERAll': 0.0, '65EGA': 0.0, '65EGB': 0.0, '65EGC': 0.0, '65EGD': 0.0, '65EGF': 0.0, '65EGG': 0.0, '65EGCRF01_AE': 0.1, '65EGCRF02_AG': 0.0, '65EGOther': 0.1, '65EGAll': 0.0, '65EVA': 0.0, '65EVB': 0.0, '65EVC': 0.0, '65EVD': 0.0, '65EVF': 0.1, '65EVG': 0.0, '65EVCRF01_AE': 0.0, '65EVCRF02_AG': 0.0, '65EVOther': 0.0, '65EVAll': 0.0, '65EHA': 0.0, '65EHB': 0.0, '65EHC': 0.0, '65EHD': 0.0, '65EHF': 0.1, '65EHG': 0.0, '65EHCRF01_AE': 0.0, '65EHCRF02_AG': 0.0, '65EHOther': 0.0, '65EHAll': 0.0, '65ETA': 0.0, '65ETB': 0.0, '65ETC': 0.0, '65ETD': 0.0, '65ETF': 0.0, '65ETG': 0.0, '65ETCRF01_AE': 0.0, '65ETCRF02_AG': 0.0, '65ETOther': 0.0, '65ETAll': 0.0, '65EAA': 0.0, '65EAB': 0.0, '65EAC': 0.0, '65EAD': 0.0, '65EAF': 0.0, '65EAG': 0.0, '65EACRF01_AE': 0.0, '65EACRF02_AG': 0.0, '65EAOther': 0.0, '65EAAll': 0.0, '65EIA': 0.0, '65EIB': 0.0, '65EIC': 0.0, '65EID': 0.0, '65EIF': 0.0, '65EIG': 0.0, '65EICRF01_AE': 0.0, '65EICRF02_AG': 0.0, '65EIOther': 0.0, '65EIAll': 0.0, '66IFA': 1.5, '66IFB': 1.4, '66IFC': 0.0, '66IFD': 1.3, '66IFF': 0.9, '66IFG': 6.8, '66IFCRF01_AE': 1.5, '66IFCRF02_AG': 8.8, '66IFOther': 0.5, '66IFAll': 1.3, '66IVA': 0.0, '66IVB': 1.4, '66IVC': 0.3, '66IVD': 0.0, '66IVF': 0.1, '66IVG': 0.8, '66IVCRF01_AE': 0.0, '66IVCRF02_AG': 0.4, '66IVOther': 0.7, '66IVAll': 1.0, '66ILA': 0.0, '66ILB': 0.3, '66ILC': 0.1, '66ILD': 0.0, '66ILF': 0.1, '66ILG': 0.5, '66ILCRF01_AE': 0.0, '66ILCRF02_AG': 0.0, '66ILOther': 0.2, '66ILAll': 0.3, '66IMA': 0.0, '66IMB': 0.0, '66IMC': 0.0, '66IMD': 0.0, '66IMF': 0.1, '66IMG': 0.0, '66IMCRF01_AE': 0.0, '66IMCRF02_AG': 0.0, '66IMOther': 0.1, '66IMAll': 0.0, '66INA': 0.3, '66INB': 0.0, '66INC': 0.0, '66IND': 0.0, '66INF': 0.0, '66ING': 0.0, '66INCRF01_AE': 0.0, '66INCRF02_AG': 0.0, '66INOther': 0.0, '66INAll': 0.0, '66ITA': 0.0, '66ITB': 0.0, '66ITC': 0.0, '66ITD': 0.0, '66ITF': 0.0, '66ITG': 0.0, '66ITCRF01_AE': 0.0, '66ITCRF02_AG': 0.0, '66ITOther': 0.1, '66ITAll': 0.0, '66ISA': 0.0, '66ISB': 0.0, '66ISC': 0.0, '66ISD': 0.0, '66ISF': 0.0, '66ISG': 0.0, '66ISCRF01_AE': 0.0, '66ISCRF02_AG': 0.0, '66ISOther': 0.0, '66ISAll': 0.0, '66IEA': 0.0, '66IEB': 0.0, '66IEC': 0.0, '66IED': 0.0, '66IEF': 0.0, '66IEG': 0.0, '66IECRF01_AE': 0.0, '66IECRF02_AG': 0.0, '66IEOther': 0.0, '66IEAll': 0.0, '66ICA': 0.0, '66ICB': 0.0, '66ICC': 0.0, '66ICD': 0.0, '66ICF': 0.0, '66ICG': 0.0, '66ICCRF01_AE': 0.0, '66ICCRF02_AG': 0.0, '66ICOther': 0.0, '66ICAll': 0.0, '67CEA': 0.0, '67CEB': 0.2, '67CEC': 0.1, '67CED': 0.4, '67CEF': 0.6, '67CEG': 11.4, '67CECRF01_AE': 0.0, '67CECRF02_AG': 0.0, '67CEOther': 1.6, '67CEAll': 0.5, '67CYA': 0.9, '67CYB': 1.0, '67CYC': 0.3, '67CYD': 1.7, '67CYF': 0.3, '67CYG': 2.3, '67CYCRF01_AE': 0.6, '67CYCRF02_AG': 2.7, '67CYOther': 0.5, '67CYAll': 0.8, '67CSA': 0.0, '67CSB': 0.8, '67CSC': 0.3, '67CSD': 0.8, '67CSF': 0.5, '67CSG': 2.3, '67CSCRF01_AE': 0.1, '67CSCRF02_AG': 0.4, '67CSOther': 0.3, '67CSAll': 0.7, '67CFA': 0.3, '67CFB': 1.2, '67CFC': 0.3, '67CFD': 0.4, '67CFF': 0.1, '67CFG': 1.0, '67CFCRF01_AE': 0.1, '67CFCRF02_AG': 0.0, '67CFOther': 0.3, '67CFAll': 0.9, '67CDA': 0.0, '67CDB': 0.3, '67CDC': 0.0, '67CDD': 0.0, '67CDF': 0.6, '67CDG': 0.3, '67CDCRF01_AE': 0.1, '67CDCRF02_AG': 0.0, '67CDOther': 0.7, '67CDAll': 0.3, '67CGA': 0.0, '67CGB': 0.1, '67CGC': 0.0, '67CGD': 0.0, '67CGF': 0.0, '67CGG': 0.5, '67CGCRF01_AE': 0.3, '67CGCRF02_AG': 0.0, '67CGOther': 0.1, '67CGAll': 0.1, '67CHA': 0.0, '67CHB': 0.0, '67CHC': 0.1, '67CHD': 0.0, '67CHF': 0.1, '67CHG': 0.5, '67CHCRF01_AE': 0.0, '67CHCRF02_AG': 0.8, '67CHOther': 0.1, '67CHAll': 0.1, '67CWA': 0.0, '67CWB': 0.2, '67CWC': 0.0, '67CWD': 0.0, '67CWF': 0.2, '67CWG': 0.8, '67CWCRF01_AE': 0.1, '67CWCRF02_AG': 0.0, '67CWOther': 0.1, '67CWAll': 0.2, '67CNA': 0.0, '67CNB': 0.0, '67CNC': 0.0, '67CND': 0.0, '67CNF': 0.0, '67CNG': 0.8, '67CNCRF01_AE': 0.0, '67CNCRF02_AG': 0.0, '67CNOther': 0.0, '67CNAll': 0.0, '67CQA': 0.0, '67CQB': 0.0, '67CQC': 0.0, '67CQD': 0.0, '67CQF': 0.1, '67CQG': 0.0, '67CQCRF01_AE': 0.0, '67CQCRF02_AG': 0.0, '67CQOther': 0.2, '67CQAll': 0.0, '67CRA': 0.0, '67CRB': 0.0, '67CRC': 0.0, '67CRD': 0.0, '67CRF': 0.0, '67CRG': 0.3, '67CRCRF01_AE': 0.0, '67CRCRF02_AG': 0.0, '67CROther': 0.1, '67CRAll': 0.0, '67CLA': 0.0, '67CLB': 0.2, '67CLC': 0.0, '67CLD': 0.0, '67CLF': 0.0, '67CLG': 0.0, '67CLCRF01_AE': 0.0, '67CLCRF02_AG': 0.0, '67CLOther': 0.0, '67CLAll': 0.1, '67CMA': 0.0, '67CMB': 0.0, '67CMC': 0.0, '67CMD': 0.0, '67CMF': 0.0, '67CMG': 0.0, '67CMCRF01_AE': 0.0, '67CMCRF02_AG': 0.0, '67CMOther': 0.0, '67CMAll': 0.0, '67CAA': 0.0, '67CAB': 0.0, '67CAC': 0.0, '67CAD': 0.0, '67CAF': 0.0, '67CAG': 0.0, '67CACRF01_AE': 0.0, '67CACRF02_AG': 0.0, '67CAOther': 0.0, '67CAAll': 0.0, '67CKA': 0.0, '67CKB': 0.0, '67CKC': 0.0, '67CKD': 0.0, '67CKF': 0.0, '67CKG': 0.0, '67CKCRF01_AE': 0.0, '67CKCRF02_AG': 0.0, '67CKOther': 0.0, '67CKAll': 0.0, '67CVA': 0.0, '67CVB': 0.0, '67CVC': 0.0, '67CVD': 0.0, '67CVF': 0.0, '67CVG': 0.0, '67CVCRF01_AE': 0.0, '67CVCRF02_AG': 0.0, '67CVOther': 0.0, '67CVAll': 0.0, '67CIA': 0.0, '67CIB': 0.0, '67CIC': 0.0, '67CID': 0.0, '67CIF': 0.0, '67CIG': 0.0, '67CICRF01_AE': 0.0, '67CICRF02_AG': 0.0, '67CIOther': 0.0, '67CIAll': 0.0, '68GEA': 2.1, '68GEB': 0.4, '68GEC': 0.8, '68GED': 0.0, '68GEF': 1.6, '68GEG': 1.3, '68GECRF01_AE': 0.3, '68GECRF02_AG': 1.9, '68GEOther': 1.1, '68GEAll': 0.6, '68GDA': 0.0, '68GDB': 0.1, '68GDC': 0.0, '68GDD': 0.0, '68GDF': 0.1, '68GDG': 0.5, '68GDCRF01_AE': 0.0, '68GDCRF02_AG': 0.0, '68GDOther': 0.2, '68GDAll': 0.1, '68GRA': 0.0, '68GRB': 0.0, '68GRC': 0.1, '68GRD': 0.0, '68GRF': 0.0, '68GRG': 0.0, '68GRCRF01_AE': 0.0, '68GRCRF02_AG': 0.0, '68GROther': 0.0, '68GRAll': 0.0, '68GAA': 0.0, '68GAB': 0.0, '68GAC': 0.0, '68GAD': 0.0, '68GAF': 0.0, '68GAG': 0.0, '68GACRF01_AE': 0.0, '68GACRF02_AG': 0.0, '68GAOther': 0.0, '68GAAll': 0.0, '68GKA': 0.0, '68GKB': 0.0, '68GKC': 0.0, '68GKD': 0.0, '68GKF': 0.0, '68GKG': 0.0, '68GKCRF01_AE': 0.0, '68GKCRF02_AG': 0.0, '68GKOther': 0.0, '68GKAll': 0.0, '68GVA': 0.0, '68GVB': 0.0, '68GVC': 0.0, '68GVD': 0.0, '68GVF': 0.0, '68GVG': 0.0, '68GVCRF01_AE': 0.0, '68GVCRF02_AG': 0.0, '68GVOther': 0.0, '68GVAll': 0.0, '68GQA': 0.0, '68GQB': 0.0, '68GQC': 0.0, '68GQD': 0.0, '68GQF': 0.0, '68GQG': 0.0, '68GQCRF01_AE': 0.0, '68GQCRF02_AG': 0.0, '68GQOther': 0.0, '68GQAll': 0.0, '68GCA': 0.0, '68GCB': 0.0, '68GCC': 0.0, '68GCD': 0.0, '68GCF': 0.0, '68GCG': 0.0, '68GCCRF01_AE': 0.0, '68GCCRF02_AG': 0.0, '68GCOther': 0.0, '68GCAll': 0.0, '69HKA': 97.3, '69HKB': 3.1, '69HKC': 98.5, '69HKD': 4.6, '69HKF': 0.1, '69HKG': 96.4, '69HKCRF01_AE': 97.3, '69HKCRF02_AG': 95.4, '69HKOther': 9.6, '69HKAll': 20.7, '69HQA': 1.2, '69HQB': 4.7, '69HQC': 1.0, '69HQD': 20.9, '69HQF': 1.0, '69HQG': 0.3, '69HQCRF01_AE': 0.4, '69HQCRF02_AG': 1.9, '69HQOther': 1.3, '69HQAll': 3.7, '69HYA': 0.0, '69HYB': 2.4, '69HYC': 0.0, '69HYD': 7.1, '69HYF': 4.2, '69HYG': 0.0, '69HYCRF01_AE': 0.0, '69HYCRF02_AG': 0.0, '69HYOther': 3.5, '69HYAll': 2.2, '69HRA': 0.9, '69HRB': 1.6, '69HRC': 0.4, '69HRD': 2.1, '69HRF': 1.8, '69HRG': 2.3, '69HRCRF01_AE': 0.4, '69HRCRF02_AG': 1.2, '69HROther': 0.7, '69HRAll': 1.4, '69HNA': 0.0, '69HNB': 0.5, '69HNC': 0.0, '69HND': 0.0, '69HNF': 0.1, '69HNG': 0.0, '69HNCRF01_AE': 0.1, '69HNCRF02_AG': 0.0, '69HNOther': 0.2, '69HNAll': 0.4, '69HCA': 0.0, '69HCB': 0.0, '69HCC': 0.0, '69HCD': 0.4, '69HCF': 0.2, '69HCG': 0.0, '69HCCRF01_AE': 0.0, '69HCCRF02_AG': 0.0, '69HCOther': 0.2, '69HCAll': 0.1, '69HSA': 0.0, '69HSB': 0.0, '69HSC': 0.0, '69HSD': 0.0, '69HSF': 0.1, '69HSG': 0.0, '69HSCRF01_AE': 0.0, '69HSCRF02_AG': 0.4, '69HSOther': 0.0, '69HSAll': 0.0, '69HEA': 0.0, '69HEB': 0.0, '69HEC': 0.0, '69HED': 0.0, '69HEF': 0.0, '69HEG': 0.5, '69HECRF01_AE': 0.1, '69HECRF02_AG': 0.0, '69HEOther': 0.0, '69HEAll': 0.0, '69HIA': 0.0, '69HIB': 0.1, '69HIC': 0.0, '69HID': 0.0, '69HIF': 0.0, '69HIG': 0.0, '69HICRF01_AE': 0.0, '69HICRF02_AG': 0.0, '69HIOther': 0.1, '69HIAll': 0.0, '69HTA': 0.0, '69HTB': 0.0, '69HTC': 0.0, '69HTD': 0.0, '69HTF': 0.0, '69HTG': 0.0, '69HTCRF01_AE': 0.0, '69HTCRF02_AG': 0.0, '69HTOther': 0.0, '69HTAll': 0.0, '69HLA': 0.0, '69HLB': 0.0, '69HLC': 0.0, '69HLD': 0.0, '69HLF': 0.0, '69HLG': 0.0, '69HLCRF01_AE': 0.0, '69HLCRF02_AG': 0.0, '69HLOther': 0.1, '69HLAll': 0.0, '69HAA': 0.0, '69HAB': 0.0, '69HAC': 0.0, '69HAD': 0.0, '69HAF': 0.0, '69HAG': 0.0, '69HACRF01_AE': 0.0, '69HACRF02_AG': 0.0, '69HAOther': 0.0, '69HAAll': 0.0, '69HFA': 0.0, '69HFB': 0.0, '69HFC': 0.0, '69HFD': 0.0, '69HFF': 0.0, '69HFG': 0.0, '69HFCRF01_AE': 0.0, '69HFCRF02_AG': 0.0, '69HFOther': 0.0, '69HFAll': 0.0, '69HVA': 0.0, '69HVB': 0.0, '69HVC': 0.0, '69HVD': 0.0, '69HVF': 0.0, '69HVG': 0.0, '69HVCRF01_AE': 0.0, '69HVCRF02_AG': 0.0, '69HVOther': 0.0, '69HVAll': 0.0, '69HMA': 0.0, '69HMB': 0.0, '69HMC': 0.0, '69HMD': 0.0, '69HMF': 0.0, '69HMG': 0.0, '69HMCRF01_AE': 0.0, '69HMCRF02_AG': 0.0, '69HMOther': 0.0, '69HMAll': 0.0, '69H#A': 0.0, '69H#B': 0.0, '69H#C': 0.0, '69H#D': 0.0, '69H#F': 0.0, '69H#G': 0.0, '69H#CRF01_AE': 0.0, '69H#CRF02_AG': 0.0, '69H#Other': 0.0, '69H#All': 0.0, '69HPA': 0.0, '69HPB': 0.0, '69HPC': 0.0, '69HPD': 0.0, '69HPF': 0.0, '69HPG': 0.0, '69HPCRF01_AE': 0.0, '69HPCRF02_AG': 0.0, '69HPOther': 0.0, '69HPAll': 0.0, '69HDA': 0.0, '69HDB': 0.0, '69HDC': 0.0, '69HDD': 0.0, '69HDF': 0.0, '69HDG': 0.0, '69HDCRF01_AE': 0.0, '69HDCRF02_AG': 0.0, '69HDOther': 0.0, '69HDAll': 0.0, '70KRA': 12.8, '70KRB': 3.7, '70KRC': 3.3, '70KRD': 5.4, '70KRF': 4.1, '70KRG': 17.4, '70KRCRF01_AE': 15.9, '70KRCRF02_AG': 34.5, '70KROther': 3.5, '70KRAll': 4.6, '70KEA': 0.0, '70KEB': 0.9, '70KEC': 0.1, '70KED': 1.3, '70KEF': 0.4, '70KEG': 0.3, '70KECRF01_AE': 0.0, '70KECRF02_AG': 1.1, '70KEOther': 0.7, '70KEAll': 0.7, '70KTA': 0.0, '70KTB': 0.8, '70KTC': 0.1, '70KTD': 1.7, '70KTF': 0.3, '70KTG': 0.3, '70KTCRF01_AE': 0.1, '70KTCRF02_AG': 0.0, '70KTOther': 0.3, '70KTAll': 0.6, '70KQA': 0.3, '70KQB': 0.3, '70KQC': 0.5, '70KQD': 1.3, '70KQF': 0.3, '70KQG': 0.3, '70KQCRF01_AE': 0.0, '70KQCRF02_AG': 0.4, '70KQOther': 0.3, '70KQAll': 0.3, '70KNA': 0.0, '70KNB': 0.1, '70KNC': 0.0, '70KND': 0.0, '70KNF': 0.1, '70KNG': 0.0, '70KNCRF01_AE': 0.1, '70KNCRF02_AG': 0.4, '70KNOther': 0.1, '70KNAll': 0.1, '70KFA': 0.0, '70KFB': 0.0, '70KFC': 0.0, '70KFD': 0.4, '70KFF': 0.0, '70KFG': 0.0, '70KFCRF01_AE': 0.0, '70KFCRF02_AG': 0.0, '70KFOther': 0.0, '70KFAll': 0.0, '70KIA': 0.0, '70KIB': 0.1, '70KIC': 0.0, '70KID': 0.0, '70KIF': 0.0, '70KIG': 0.0, '70KICRF01_AE': 0.0, '70KICRF02_AG': 0.0, '70KIOther': 0.0, '70KIAll': 0.0, '70KMA': 0.0, '70KMB': 0.0, '70KMC': 0.0, '70KMD': 0.0, '70KMF': 0.0, '70KMG': 0.0, '70KMCRF01_AE': 0.0, '70KMCRF02_AG': 0.0, '70KMOther': 0.0, '70KMAll': 0.0, '70KGA': 0.0, '70KGB': 0.0, '70KGC': 0.0, '70KGD': 0.0, '70KGF': 0.0, '70KGG': 0.0, '70KGCRF01_AE': 0.0, '70KGCRF02_AG': 0.0, '70KGOther': 0.0, '70KGAll': 0.0, '70KSA': 0.0, '70KSB': 0.0, '70KSC': 0.0, '70KSD': 0.0, '70KSF': 0.0, '70KSG': 0.0, '70KSCRF01_AE': 0.0, '70KSCRF02_AG': 0.0, '70KSOther': 0.0, '70KSAll': 0.0, '70KYA': 0.0, '70KYB': 0.0, '70KYC': 0.0, '70KYD': 0.0, '70KYF': 0.0, '70KYG': 0.0, '70KYCRF01_AE': 0.0, '70KYCRF02_AG': 0.0, '70KYOther': 0.0, '70KYAll': 0.0, '70KVA': 0.0, '70KVB': 0.0, '70KVC': 0.0, '70KVD': 0.0, '70KVF': 0.0, '70KVG': 0.0, '70KVCRF01_AE': 0.0, '70KVCRF02_AG': 0.0, '70KVOther': 0.0, '70KVAll': 0.0, '70KLA': 0.0, '70KLB': 0.0, '70KLC': 0.0, '70KLD': 0.0, '70KLF': 0.0, '70KLG': 0.0, '70KLCRF01_AE': 0.0, '70KLCRF02_AG': 0.0, '70KLOther': 0.0, '70KLAll': 0.0, '70KAA': 0.0, '70KAB': 0.0, '70KAC': 0.0, '70KAD': 0.0, '70KAF': 0.0, '70KAG': 0.0, '70KACRF01_AE': 0.0, '70KACRF02_AG': 0.0, '70KAOther': 0.0, '70KAAll': 0.0, '70KPA': 0.0, '70KPB': 0.0, '70KPC': 0.0, '70KPD': 0.0, '70KPF': 0.0, '70KPG': 0.0, '70KPCRF01_AE': 0.0, '70KPCRF02_AG': 0.0, '70KPOther': 0.0, '70KPAll': 0.0, '70KHA': 0.0, '70KHB': 0.0, '70KHC': 0.0, '70KHD': 0.0, '70KHF': 0.0, '70KHG': 0.0, '70KHCRF01_AE': 0.0, '70KHCRF02_AG': 0.0, '70KHOther': 0.0, '70KHAll': 0.0, '71AVA': 5.6, '71AVB': 28.3, '71AVC': 7.6, '71AVD': 18.1, '71AVF': 14.5, '71AVG': 18.0, '71AVCRF01_AE': 3.6, '71AVCRF02_AG': 8.0, '71AVOther': 18.5, '71AVAll': 22.7, '71ATA': 0.9, '71ATB': 10.6, '71ATC': 1.6, '71ATD': 5.9, '71ATF': 4.5, '71ATG': 6.0, '71ATCRF01_AE': 1.9, '71ATCRF02_AG': 3.8, '71ATOther': 4.8, '71ATAll': 8.2, '71AIA': 0.6, '71AIB': 3.2, '71AIC': 0.5, '71AID': 2.5, '71AIF': 1.5, '71AIG': 1.0, '71AICRF01_AE': 0.3, '71AICRF02_AG': 0.4, '71AIOther': 0.8, '71AIAll': 2.6, '71ALA': 0.0, '71ALB': 0.5, '71ALC': 0.0, '71ALD': 0.0, '71ALF': 0.2, '71ALG': 0.3, '71ALCRF01_AE': 0.3, '71ALCRF02_AG': 0.0, '71ALOther': 0.3, '71ALAll': 0.4, '71AGA': 0.0, '71AGB': 0.0, '71AGC': 0.0, '71AGD': 0.0, '71AGF': 0.0, '71AGG': 0.0, '71AGCRF01_AE': 0.0, '71AGCRF02_AG': 0.4, '71AGOther': 0.0, '71AGAll': 0.0, '71AMA': 0.0, '71AMB': 0.0, '71AMC': 0.0, '71AMD': 0.4, '71AMF': 0.0, '71AMG': 0.0, '71AMCRF01_AE': 0.0, '71AMCRF02_AG': 0.0, '71AMOther': 0.0, '71AMAll': 0.0, '71ADA': 0.0, '71ADB': 0.0, '71ADC': 0.0, '71ADD': 0.0, '71ADF': 0.0, '71ADG': 0.3, '71ADCRF01_AE': 0.0, '71ADCRF02_AG': 0.0, '71ADOther': 0.0, '71ADAll': 0.0, '71ASA': 0.0, '71ASB': 0.0, '71ASC': 0.0, '71ASD': 0.0, '71ASF': 0.0, '71ASG': 0.0, '71ASCRF01_AE': 0.0, '71ASCRF02_AG': 0.0, '71ASOther': 0.1, '71ASAll': 0.0, '71ANA': 0.0, '71ANB': 0.0, '71ANC': 0.0, '71AND': 0.0, '71ANF': 0.0, '71ANG': 0.0, '71ANCRF01_AE': 0.0, '71ANCRF02_AG': 0.0, '71ANOther': 0.0, '71ANAll': 0.0, '71APA': 0.0, '71APB': 0.0, '71APC': 0.0, '71APD': 0.0, '71APF': 0.0, '71APG': 0.0, '71APCRF01_AE': 0.0, '71APCRF02_AG': 0.0, '71APOther': 0.0, '71APAll': 0.0, '71AFA': 0.0, '71AFB': 0.0, '71AFC': 0.0, '71AFD': 0.0, '71AFF': 0.0, '71AFG': 0.0, '71AFCRF01_AE': 0.0, '71AFCRF02_AG': 0.0, '71AFOther': 0.0, '71AFAll': 0.0, '71ACA': 0.0, '71ACB': 0.0, '71ACC': 0.0, '71ACD': 0.0, '71ACF': 0.0, '71ACG': 0.0, '71ACCRF01_AE': 0.0, '71ACCRF02_AG': 0.0, '71ACOther': 0.0, '71ACAll': 0.0, '71A#A': 0.0, '71A#B': 0.0, '71A#C': 0.0, '71A#D': 0.0, '71A#F': 0.0, '71A#G': 0.0, '71A#CRF01_AE': 0.0, '71A#CRF02_AG': 0.0, '71A#Other': 0.0, '71A#All': 0.0, '72ITA': 1.2, '72ITB': 6.2, '72ITC': 0.9, '72ITD': 2.9, '72ITF': 31.3, '72ITG': 2.1, '72ITCRF01_AE': 1.5, '72ITCRF02_AG': 1.1, '72ITOther': 24.0, '72ITAll': 8.1, '72IVA': 2.7, '72IVB': 12.1, '72IVC': 2.9, '72IVD': 15.5, '72IVF': 8.9, '72IVG': 2.3, '72IVCRF01_AE': 6.1, '72IVCRF02_AG': 6.5, '72IVOther': 4.0, '72IVAll': 9.8, '72IMA': 0.0, '72IMB': 1.9, '72IMC': 0.8, '72IMD': 1.7, '72IMF': 0.7, '72IMG': 0.5, '72IMCRF01_AE': 0.1, '72IMCRF02_AG': 2.7, '72IMOther': 1.2, '72IMAll': 1.5, '72IRA': 0.6, '72IRB': 0.9, '72IRC': 0.2, '72IRD': 0.0, '72IRF': 3.5, '72IRG': 0.3, '72IRCRF01_AE': 0.6, '72IRCRF02_AG': 0.4, '72IROther': 3.7, '72IRAll': 1.1, '72IEA': 0.0, '72IEB': 1.1, '72IEC': 0.0, '72IED': 0.4, '72IEF': 0.9, '72IEG': 0.3, '72IECRF01_AE': 0.1, '72IECRF02_AG': 0.0, '72IEOther': 0.8, '72IEAll': 0.9, '72ILA': 0.0, '72ILB': 2.8, '72ILC': 0.5, '72ILD': 0.0, '72ILF': 0.3, '72ILG': 0.0, '72ILCRF01_AE': 0.0, '72ILCRF02_AG': 0.0, '72ILOther': 0.2, '72ILAll': 2.0, '72IKA': 0.3, '72IKB': 0.5, '72IKC': 0.1, '72IKD': 0.0, '72IKF': 0.2, '72IKG': 0.0, '72IKCRF01_AE': 0.1, '72IKCRF02_AG': 0.0, '72IKOther': 0.8, '72IKAll': 0.4, '72ICA': 0.0, '72ICB': 0.0, '72ICC': 0.0, '72ICD': 0.4, '72ICF': 0.0, '72ICG': 0.0, '72ICCRF01_AE': 0.0, '72ICCRF02_AG': 0.4, '72ICOther': 0.0, '72ICAll': 0.0, '72ISA': 0.0, '72ISB': 0.0, '72ISC': 0.0, '72ISD': 0.4, '72ISF': 0.1, '72ISG': 0.0, '72ISCRF01_AE': 0.0, '72ISCRF02_AG': 0.0, '72ISOther': 0.1, '72ISAll': 0.0, '72IQA': 0.0, '72IQB': 0.0, '72IQC': 0.0, '72IQD': 0.0, '72IQF': 0.1, '72IQG': 0.0, '72IQCRF01_AE': 0.0, '72IQCRF02_AG': 0.0, '72IQOther': 0.1, '72IQAll': 0.0, '72IFA': 0.0, '72IFB': 0.1, '72IFC': 0.0, '72IFD': 0.0, '72IFF': 0.0, '72IFG': 0.0, '72IFCRF01_AE': 0.0, '72IFCRF02_AG': 0.0, '72IFOther': 0.0, '72IFAll': 0.0, '72IAA': 0.0, '72IAB': 0.0, '72IAC': 0.0, '72IAD': 0.0, '72IAF': 0.0, '72IAG': 0.0, '72IACRF01_AE': 0.0, '72IACRF02_AG': 0.0, '72IAOther': 0.0, '72IAAll': 0.0, '72IPA': 0.0, '72IPB': 0.0, '72IPC': 0.0, '72IPD': 0.0, '72IPF': 0.1, '72IPG': 0.0, '72IPCRF01_AE': 0.0, '72IPCRF02_AG': 0.0, '72IPOther': 0.0, '72IPAll': 0.0, '72INA': 0.0, '72INB': 0.0, '72INC': 0.0, '72IND': 0.0, '72INF': 0.0, '72ING': 0.0, '72INCRF01_AE': 0.0, '72INCRF02_AG': 0.0, '72INOther': 0.0, '72INAll': 0.0, '72IYA': 0.0, '72IYB': 0.0, '72IYC': 0.0, '72IYD': 0.0, '72IYF': 0.0, '72IYG': 0.0, '72IYCRF01_AE': 0.0, '72IYCRF02_AG': 0.0, '72IYOther': 0.0, '72IYAll': 0.0, '72IWA': 0.0, '72IWB': 0.0, '72IWC': 0.0, '72IWD': 0.0, '72IWF': 0.0, '72IWG': 0.0, '72IWCRF01_AE': 0.0, '72IWCRF02_AG': 0.0, '72IWOther': 0.0, '72IWAll': 0.0, '72IDA': 0.0, '72IDB': 0.0, '72IDC': 0.0, '72IDD': 0.0, '72IDF': 0.0, '72IDG': 0.0, '72IDCRF01_AE': 0.0, '72IDCRF02_AG': 0.0, '72IDOther': 0.0, '72IDAll': 0.0, '73GSA': 1.8, '73GSB': 8.8, '73GSC': 1.0, '73GSD': 5.0, '73GSF': 2.0, '73GSG': 0.0, '73GSCRF01_AE': 0.0, '73GSCRF02_AG': 1.5, '73GSOther': 3.8, '73GSAll': 6.5, '73GTA': 0.0, '73GTB': 2.7, '73GTC': 0.2, '73GTD': 1.3, '73GTF': 0.2, '73GTG': 0.0, '73GTCRF01_AE': 0.0, '73GTCRF02_AG': 0.0, '73GTOther': 0.3, '73GTAll': 1.9, '73GCA': 0.0, '73GCB': 1.2, '73GCC': 0.1, '73GCD': 0.4, '73GCF': 0.4, '73GCG': 0.3, '73GCCRF01_AE': 0.0, '73GCCRF02_AG': 0.0, '73GCOther': 0.9, '73GCAll': 0.9, '73GDA': 0.3, '73GDB': 0.1, '73GDC': 0.0, '73GDD': 0.4, '73GDF': 0.0, '73GDG': 0.3, '73GDCRF01_AE': 0.0, '73GDCRF02_AG': 0.4, '73GDOther': 0.1, '73GDAll': 0.1, '73GAA': 0.0, '73GAB': 0.5, '73GAC': 0.0, '73GAD': 0.0, '73GAF': 0.1, '73GAG': 0.3, '73GACRF01_AE': 0.0, '73GACRF02_AG': 0.4, '73GAOther': 0.1, '73GAAll': 0.4, '73GVA': 0.0, '73GVB': 0.1, '73GVC': 0.0, '73GVD': 0.0, '73GVF': 0.0, '73GVG': 0.0, '73GVCRF01_AE': 0.0, '73GVCRF02_AG': 0.0, '73GVOther': 0.1, '73GVAll': 0.1, '73GEA': 0.0, '73GEB': 0.0, '73GEC': 0.0, '73GED': 0.0, '73GEF': 0.0, '73GEG': 0.3, '73GECRF01_AE': 0.0, '73GECRF02_AG': 0.0, '73GEOther': 0.0, '73GEAll': 0.0, '73GNA': 0.0, '73GNB': 0.0, '73GNC': 0.0, '73GND': 0.0, '73GNF': 0.1, '73GNG': 0.0, '73GNCRF01_AE': 0.0, '73GNCRF02_AG': 0.0, '73GNOther': 0.0, '73GNAll': 0.0, '73GRA': 0.0, '73GRB': 0.0, '73GRC': 0.0, '73GRD': 0.0, '73GRF': 0.0, '73GRG': 0.0, '73GRCRF01_AE': 0.0, '73GRCRF02_AG': 0.0, '73GROther': 0.1, '73GRAll': 0.0, '73GKA': 0.0, '73GKB': 0.0, '73GKC': 0.0, '73GKD': 0.0, '73GKF': 0.0, '73GKG': 0.0, '73GKCRF01_AE': 0.0, '73GKCRF02_AG': 0.0, '73GKOther': 0.0, '73GKAll': 0.0, '73GIA': 0.0, '73GIB': 0.0, '73GIC': 0.0, '73GID': 0.0, '73GIF': 0.0, '73GIG': 0.0, '73GICRF01_AE': 0.0, '73GICRF02_AG': 0.0, '73GIOther': 0.0, '73GIAll': 0.0, '74TSA': 6.8, '74TSB': 5.2, '74TSC': 16.5, '74TSD': 7.9, '74TSF': 13.8, '74TSG': 15.8, '74TSCRF01_AE': 7.1, '74TSCRF02_AG': 9.6, '74TSOther': 11.0, '74TSAll': 7.8, '74TAA': 2.7, '74TAB': 1.2, '74TAC': 2.6, '74TAD': 0.8, '74TAF': 2.9, '74TAG': 1.3, '74TACRF01_AE': 2.4, '74TACRF02_AG': 2.3, '74TAOther': 1.2, '74TAAll': 1.5, '74TPA': 1.8, '74TPB': 1.8, '74TPC': 2.1, '74TPD': 1.7, '74TPF': 3.4, '74TPG': 0.8, '74TPCRF01_AE': 0.7, '74TPCRF02_AG': 0.4, '74TPOther': 0.9, '74TPAll': 1.8, '74TKA': 0.0, '74TKB': 0.2, '74TKC': 0.2, '74TKD': 0.0, '74TKF': 0.0, '74TKG': 0.3, '74TKCRF01_AE': 0.0, '74TKCRF02_AG': 0.0, '74TKOther': 0.1, '74TKAll': 0.2, '74TEA': 0.0, '74TEB': 0.1, '74TEC': 0.2, '74TED': 0.0, '74TEF': 0.1, '74TEG': 0.0, '74TECRF01_AE': 0.0, '74TECRF02_AG': 0.4, '74TEOther': 0.1, '74TEAll': 0.1, '74TIA': 0.0, '74TIB': 0.0, '74TIC': 0.0, '74TID': 0.0, '74TIF': 0.0, '74TIG': 0.0, '74TICRF01_AE': 0.0, '74TICRF02_AG': 0.0, '74TIOther': 0.1, '74TIAll': 0.0, '74TMA': 0.0, '74TMB': 0.0, '74TMC': 0.0, '74TMD': 0.0, '74TMF': 0.0, '74TMG': 0.0, '74TMCRF01_AE': 0.0, '74TMCRF02_AG': 0.0, '74TMOther': 0.1, '74TMAll': 0.0, '74TRA': 0.0, '74TRB': 0.0, '74TRC': 0.1, '74TRD': 0.0, '74TRF': 0.0, '74TRG': 0.0, '74TRCRF01_AE': 0.0, '74TRCRF02_AG': 0.0, '74TROther': 0.0, '74TRAll': 0.0, '74TNA': 0.0, '74TNB': 0.0, '74TNC': 0.0, '74TND': 0.0, '74TNF': 0.0, '74TNG': 0.0, '74TNCRF01_AE': 0.0, '74TNCRF02_AG': 0.0, '74TNOther': 0.0, '74TNAll': 0.0, '74TVA': 0.0, '74TVB': 0.0, '74TVC': 0.0, '74TVD': 0.0, '74TVF': 0.0, '74TVG': 0.0, '74TVCRF01_AE': 0.0, '74TVCRF02_AG': 0.0, '74TVOther': 0.0, '74TVAll': 0.0, '74TQA': 0.0, '74TQB': 0.0, '74TQC': 0.0, '74TQD': 0.0, '74TQF': 0.0, '74TQG': 0.0, '74TQCRF01_AE': 0.0, '74TQCRF02_AG': 0.0, '74TQOther': 0.0, '74TQAll': 0.0, '74TLA': 0.0, '74TLB': 0.0, '74TLC': 0.0, '74TLD': 0.0, '74TLF': 0.0, '74TLG': 0.0, '74TLCRF01_AE': 0.0, '74TLCRF02_AG': 0.0, '74TLOther': 0.0, '74TLAll': 0.0, '75VIA': 0.0, '75VIB': 1.1, '75VIC': 0.8, '75VID': 0.8, '75VIF': 0.7, '75VIG': 0.0, '75VICRF01_AE': 0.0, '75VICRF02_AG': 0.8, '75VIOther': 1.1, '75VIAll': 1.0, '75VAA': 0.0, '75VAB': 0.0, '75VAC': 0.0, '75VAD': 0.0, '75VAF': 0.1, '75VAG': 0.0, '75VACRF01_AE': 0.1, '75VACRF02_AG': 0.0, '75VAOther': 0.1, '75VAAll': 0.0, '75VMA': 0.3, '75VMB': 0.0, '75VMC': 0.0, '75VMD': 0.0, '75VMF': 0.0, '75VMG': 0.0, '75VMCRF01_AE': 0.0, '75VMCRF02_AG': 0.0, '75VMOther': 0.0, '75VMAll': 0.0, '75VGA': 0.0, '75VGB': 0.0, '75VGC': 0.0, '75VGD': 0.0, '75VGF': 0.0, '75VGG': 0.0, '75VGCRF01_AE': 0.0, '75VGCRF02_AG': 0.0, '75VGOther': 0.0, '75VGAll': 0.0, '75VFA': 0.0, '75VFB': 0.0, '75VFC': 0.0, '75VFD': 0.0, '75VFF': 0.0, '75VFG': 0.0, '75VFCRF01_AE': 0.0, '75VFCRF02_AG': 0.0, '75VFOther': 0.0, '75VFAll': 0.0, '75VEA': 0.0, '75VEB': 0.0, '75VEC': 0.0, '75VED': 0.0, '75VEF': 0.0, '75VEG': 0.0, '75VECRF01_AE': 0.0, '75VECRF02_AG': 0.0, '75VEOther': 0.0, '75VEAll': 0.0, '75VLA': 0.0, '75VLB': 0.0, '75VLC': 0.0, '75VLD': 0.0, '75VLF': 0.0, '75VLG': 0.0, '75VLCRF01_AE': 0.0, '75VLCRF02_AG': 0.0, '75VLOther': 0.0, '75VLAll': 0.0, '76LVA': 2.1, '76LVB': 3.5, '76LVC': 3.0, '76LVD': 3.8, '76LVF': 5.8, '76LVG': 2.6, '76LVCRF01_AE': 1.5, '76LVCRF02_AG': 8.8, '76LVOther': 3.3, '76LVAll': 3.5, '76LMA': 0.0, '76LMB': 0.0, '76LMC': 0.0, '76LMD': 0.0, '76LMF': 0.7, '76LMG': 0.0, '76LMCRF01_AE': 0.0, '76LMCRF02_AG': 0.0, '76LMOther': 0.1, '76LMAll': 0.0, '76LFA': 0.0, '76LFB': 0.0, '76LFC': 0.0, '76LFD': 0.0, '76LFF': 0.0, '76LFG': 0.0, '76LFCRF01_AE': 0.0, '76LFCRF02_AG': 0.0, '76LFOther': 0.0, '76LFAll': 0.0, '76LIA': 0.0, '76LIB': 0.0, '76LIC': 0.0, '76LID': 0.0, '76LIF': 0.0, '76LIG': 0.0, '76LICRF01_AE': 0.0, '76LICRF02_AG': 0.0, '76LIOther': 0.1, '76LIAll': 0.0, '76LSA': 0.0, '76LSB': 0.0, '76LSC': 0.0, '76LSD': 0.0, '76LSF': 0.0, '76LSG': 0.0, '76LSCRF01_AE': 0.0, '76LSCRF02_AG': 0.0, '76LSOther': 0.0, '76LSAll': 0.0, '76LEA': 0.0, '76LEB': 0.0, '76LEC': 0.0, '76LED': 0.0, '76LEF': 0.0, '76LEG': 0.0, '76LECRF01_AE': 0.0, '76LECRF02_AG': 0.0, '76LEOther': 0.0, '76LEAll': 0.0, '76LQA': 0.0, '76LQB': 0.0, '76LQC': 0.0, '76LQD': 0.0, '76LQF': 0.0, '76LQG': 0.0, '76LQCRF01_AE': 0.0, '76LQCRF02_AG': 0.0, '76LQOther': 0.0, '76LQAll': 0.0, '77VIA': 17.8, '77VIB': 31.7, '77VIC': 8.8, '77VID': 15.1, '77VIF': 12.1, '77VIG': 1.6, '77VICRF01_AE': 2.5, '77VICRF02_AG': 1.9, '77VIOther': 8.6, '77VIAll': 24.2, '77VLA': 0.0, '77VLB': 0.2, '77VLC': 0.0, '77VLD': 0.4, '77VLF': 0.2, '77VLG': 0.0, '77VLCRF01_AE': 0.0, '77VLCRF02_AG': 0.0, '77VLOther': 0.2, '77VLAll': 0.1, '77VTA': 0.0, '77VTB': 0.1, '77VTC': 0.0, '77VTD': 0.0, '77VTF': 0.0, '77VTG': 0.0, '77VTCRF01_AE': 0.0, '77VTCRF02_AG': 0.0, '77VTOther': 0.0, '77VTAll': 0.1, '77VEA': 0.0, '77VEB': 0.0, '77VEC': 0.0, '77VED': 0.0, '77VEF': 0.1, '77VEG': 0.0, '77VECRF01_AE': 0.1, '77VECRF02_AG': 0.0, '77VEOther': 0.0, '77VEAll': 0.0, '77VMA': 0.0, '77VMB': 0.0, '77VMC': 0.0, '77VMD': 0.0, '77VMF': 0.0, '77VMG': 0.0, '77VMCRF01_AE': 0.0, '77VMCRF02_AG': 0.0, '77VMOther': 0.0, '77VMAll': 0.0, '77VSA': 0.0, '77VSB': 0.0, '77VSC': 0.0, '77VSD': 0.0, '77VSF': 0.0, '77VSG': 0.0, '77VSCRF01_AE': 0.0, '77VSCRF02_AG': 0.0, '77VSOther': 0.0, '77VSAll': 0.0, '77VFA': 0.0, '77VFB': 0.0, '77VFC': 0.0, '77VFD': 0.0, '77VFF': 0.0, '77VFG': 0.0, '77VFCRF01_AE': 0.0, '77VFCRF02_AG': 0.0, '77VFOther': 0.0, '77VFAll': 0.0, '77VKA': 0.0, '77VKB': 0.0, '77VKC': 0.0, '77VKD': 0.0, '77VKF': 0.0, '77VKG': 0.0, '77VKCRF01_AE': 0.0, '77VKCRF02_AG': 0.0, '77VKOther': 0.0, '77VKAll': 0.0, '77VAA': 0.0, '77VAB': 0.0, '77VAC': 0.0, '77VAD': 0.0, '77VAF': 0.0, '77VAG': 0.0, '77VACRF01_AE': 0.0, '77VACRF02_AG': 0.0, '77VAOther': 0.0, '77VAAll': 0.0, '77VGA': 0.0, '77VGB': 0.0, '77VGC': 0.0, '77VGD': 0.0, '77VGF': 0.0, '77VGG': 0.0, '77VGCRF01_AE': 0.0, '77VGCRF02_AG': 0.0, '77VGOther': 0.0, '77VGAll': 0.0, '78GEA': 0.0, '78GEB': 0.0, '78GEC': 0.0, '78GED': 0.0, '78GEF': 0.0, '78GEG': 0.0, '78GECRF01_AE': 0.0, '78GECRF02_AG': 0.0, '78GEOther': 0.1, '78GEAll': 0.0, '78GRA': 0.0, '78GRB': 0.0, '78GRC': 0.0, '78GRD': 0.0, '78GRF': 0.0, '78GRG': 0.0, '78GRCRF01_AE': 0.0, '78GRCRF02_AG': 0.0, '78GROther': 0.0, '78GRAll': 0.0, '78GKA': 0.0, '78GKB': 0.0, '78GKC': 0.0, '78GKD': 0.0, '78GKF': 0.0, '78GKG': 0.0, '78GKCRF01_AE': 0.0, '78GKCRF02_AG': 0.0, '78GKOther': 0.0, '78GKAll': 0.0, '78GYA': 0.0, '78GYB': 0.0, '78GYC': 0.0, '78GYD': 0.0, '78GYF': 0.0, '78GYG': 0.0, '78GYCRF01_AE': 0.0, '78GYCRF02_AG': 0.0, '78GYOther': 0.0, '78GYAll': 0.0, '78GVA': 0.0, '78GVB': 0.0, '78GVC': 0.0, '78GVD': 0.0, '78GVF': 0.0, '78GVG': 0.0, '78GVCRF01_AE': 0.0, '78GVCRF02_AG': 0.0, '78GVOther': 0.0, '78GVAll': 0.0, '78GAA': 0.0, '78GAB': 0.0, '78GAC': 0.0, '78GAD': 0.0, '78GAF': 0.0, '78GAG': 0.0, '78GACRF01_AE': 0.0, '78GACRF02_AG': 0.0, '78GAOther': 0.0, '78GAAll': 0.0, '78GDA': 0.0, '78GDB': 0.0, '78GDC': 0.0, '78GDD': 0.0, '78GDF': 0.0, '78GDG': 0.0, '78GDCRF01_AE': 0.0, '78GDCRF02_AG': 0.0, '78GDOther': 0.0, '78GDAll': 0.0, '79PSA': 0.0, '79PSB': 0.5, '79PSC': 0.2, '79PSD': 0.0, '79PSF': 0.3, '79PSG': 0.3, '79PSCRF01_AE': 0.1, '79PSCRF02_AG': 0.8, '79PSOther': 0.4, '79PSAll': 0.4, '79PAA': 0.0, '79PAB': 0.8, '79PAC': 0.1, '79PAD': 0.4, '79PAF': 0.2, '79PAG': 0.0, '79PACRF01_AE': 0.1, '79PACRF02_AG': 0.8, '79PAOther': 0.3, '79PAAll': 0.6, '79PDA': 0.0, '79PDB': 0.4, '79PDC': 0.1, '79PDD': 0.4, '79PDF': 0.0, '79PDG': 0.0, '79PDCRF01_AE': 0.1, '79PDCRF02_AG': 0.0, '79PDOther': 0.2, '79PDAll': 0.3, '79PHA': 0.0, '79PHB': 0.1, '79PHC': 0.1, '79PHD': 0.4, '79PHF': 0.0, '79PHG': 0.0, '79PHCRF01_AE': 0.0, '79PHCRF02_AG': 0.0, '79PHOther': 0.2, '79PHAll': 0.1, '79PLA': 0.0, '79PLB': 0.0, '79PLC': 0.1, '79PLD': 0.0, '79PLF': 0.1, '79PLG': 0.0, '79PLCRF01_AE': 0.0, '79PLCRF02_AG': 0.0, '79PLOther': 0.0, '79PLAll': 0.0, '79PQA': 0.0, '79PQB': 0.1, '79PQC': 0.0, '79PQD': 0.0, '79PQF': 0.0, '79PQG': 0.0, '79PQCRF01_AE': 0.0, '79PQCRF02_AG': 0.0, '79PQOther': 0.0, '79PQAll': 0.0, '79PTA': 0.0, '79PTB': 0.0, '79PTC': 0.0, '79PTD': 0.0, '79PTF': 0.0, '79PTG': 0.0, '79PTCRF01_AE': 0.0, '79PTCRF02_AG': 0.0, '79PTOther': 0.0, '79PTAll': 0.0, '79PNA': 0.0, '79PNB': 0.1, '79PNC': 0.0, '79PND': 0.0, '79PNF': 0.0, '79PNG': 0.0, '79PNCRF01_AE': 0.0, '79PNCRF02_AG': 0.0, '79PNOther': 0.0, '79PNAll': 0.0, '79PEA': 0.0, '79PEB': 0.0, '79PEC': 0.0, '79PED': 0.0, '79PEF': 0.0, '79PEG': 0.0, '79PECRF01_AE': 0.0, '79PECRF02_AG': 0.0, '79PEOther': 0.0, '79PEAll': 0.0, '79PRA': 0.0, '79PRB': 0.0, '79PRC': 0.0, '79PRD': 0.0, '79PRF': 0.0, '79PRG': 0.0, '79PRCRF01_AE': 0.0, '79PRCRF02_AG': 0.0, '79PROther': 0.0, '79PRAll': 0.0, '79PGA': 0.0, '79PGB': 0.0, '79PGC': 0.0, '79PGD': 0.0, '79PGF': 0.0, '79PGG': 0.0, '79PGCRF01_AE': 0.0, '79PGCRF02_AG': 0.0, '79PGOther': 0.0, '79PGAll': 0.0, '80TAA': 0.3, '80TAB': 0.0, '80TAC': 0.0, '80TAD': 0.0, '80TAF': 0.0, '80TAG': 0.0, '80TACRF01_AE': 0.1, '80TACRF02_AG': 0.0, '80TAOther': 0.0, '80TAAll': 0.0, '80TVA': 0.0, '80TVB': 0.0, '80TVC': 0.0, '80TVD': 0.0, '80TVF': 0.0, '80TVG': 0.3, '80TVCRF01_AE': 0.0, '80TVCRF02_AG': 0.0, '80TVOther': 0.0, '80TVAll': 0.0, '80TQA': 0.0, '80TQB': 0.0, '80TQC': 0.0, '80TQD': 0.0, '80TQF': 0.0, '80TQG': 0.0, '80TQCRF01_AE': 0.0, '80TQCRF02_AG': 0.0, '80TQOther': 0.0, '80TQAll': 0.0, '80TIA': 0.0, '80TIB': 0.0, '80TIC': 0.0, '80TID': 0.0, '80TIF': 0.0, '80TIG': 0.0, '80TICRF01_AE': 0.0, '80TICRF02_AG': 0.0, '80TIOther': 0.1, '80TIAll': 0.0, '80TSA': 0.0, '80TSB': 0.0, '80TSC': 0.0, '80TSD': 0.0, '80TSF': 0.0, '80TSG': 0.0, '80TSCRF01_AE': 0.0, '80TSCRF02_AG': 0.0, '80TSOther': 0.0, '80TSAll': 0.0, '80TNA': 0.0, '80TNB': 0.0, '80TNC': 0.0, '80TND': 0.0, '80TNF': 0.0, '80TNG': 0.0, '80TNCRF01_AE': 0.0, '80TNCRF02_AG': 0.0, '80TNOther': 0.0, '80TNAll': 0.0, '80TKA': 0.0, '80TKB': 0.0, '80TKC': 0.0, '80TKD': 0.0, '80TKF': 0.0, '80TKG': 0.0, '80TKCRF01_AE': 0.0, '80TKCRF02_AG': 0.0, '80TKOther': 0.0, '80TKAll': 0.0, '80TPA': 0.0, '80TPB': 0.0, '80TPC': 0.0, '80TPD': 0.0, '80TPF': 0.0, '80TPG': 0.0, '80TPCRF01_AE': 0.0, '80TPCRF02_AG': 0.0, '80TPOther': 0.0, '80TPAll': 0.0, '80THA': 0.0, '80THB': 0.0, '80THC': 0.0, '80THD': 0.0, '80THF': 0.0, '80THG': 0.0, '80THCRF01_AE': 0.0, '80THCRF02_AG': 0.0, '80THOther': 0.0, '80THAll': 0.0, '80TRA': 0.0, '80TRB': 0.0, '80TRC': 0.0, '80TRD': 0.0, '80TRF': 0.0, '80TRG': 0.0, '80TRCRF01_AE': 0.0, '80TRCRF02_AG': 0.0, '80TROther': 0.0, '80TRAll': 0.0, '81PLA': 0.0, '81PLB': 0.0, '81PLC': 0.0, '81PLD': 0.0, '81PLF': 0.0, '81PLG': 0.0, '81PLCRF01_AE': 0.1, '81PLCRF02_AG': 0.4, '81PLOther': 0.0, '81PLAll': 0.0, '81PTA': 0.0, '81PTB': 0.0, '81PTC': 0.0, '81PTD': 0.0, '81PTF': 0.0, '81PTG': 0.0, '81PTCRF01_AE': 0.0, '81PTCRF02_AG': 0.0, '81PTOther': 0.0, '81PTAll': 0.0, '81PSA': 0.0, '81PSB': 0.0, '81PSC': 0.0, '81PSD': 0.0, '81PSF': 0.0, '81PSG': 0.0, '81PSCRF01_AE': 0.0, '81PSCRF02_AG': 0.0, '81PSOther': 0.0, '81PSAll': 0.0, '81PGA': 0.0, '81PGB': 0.0, '81PGC': 0.0, '81PGD': 0.0, '81PGF': 0.0, '81PGG': 0.0, '81PGCRF01_AE': 0.0, '81PGCRF02_AG': 0.0, '81PGOther': 0.1, '81PGAll': 0.0, '81PNA': 0.0, '81PNB': 0.0, '81PNC': 0.0, '81PND': 0.0, '81PNF': 0.0, '81PNG': 0.0, '81PNCRF01_AE': 0.0, '81PNCRF02_AG': 0.0, '81PNOther': 0.0, '81PNAll': 0.0, '81PAA': 0.0, '81PAB': 0.0, '81PAC': 0.0, '81PAD': 0.0, '81PAF': 0.0, '81PAG': 0.0, '81PACRF01_AE': 0.0, '81PACRF02_AG': 0.0, '81PAOther': 0.0, '81PAAll': 0.0, '81PHA': 0.0, '81PHB': 0.0, '81PHC': 0.0, '81PHD': 0.0, '81PHF': 0.0, '81PHG': 0.0, '81PHCRF01_AE': 0.0, '81PHCRF02_AG': 0.0, '81PHOther': 0.0, '81PHAll': 0.0, '81PRA': 0.0, '81PRB': 0.0, '81PRC': 0.0, '81PRD': 0.0, '81PRF': 0.0, '81PRG': 0.0, '81PRCRF01_AE': 0.0, '81PRCRF02_AG': 0.0, '81PROther': 0.0, '81PRAll': 0.0, '82VIA': 1.5, '82VIB': 2.2, '82VIC': 8.3, '82VID': 3.8, '82VIF': 3.9, '82VIG': 64.1, '82VICRF01_AE': 9.2, '82VICRF02_AG': 1.5, '82VIOther': 3.2, '82VIAll': 4.3, '82VAA': 6.5, '82VAB': 21.9, '82VAC': 11.7, '82VAD': 18.8, '82VAF': 30.2, '82VAG': 5.2, '82VACRF01_AE': 7.0, '82VACRF02_AG': 6.9, '82VAOther': 22.9, '82VAAll': 20.2, '82VTA': 0.9, '82VTB': 2.4, '82VTC': 0.2, '82VTD': 2.5, '82VTF': 1.1, '82VTG': 9.2, '82VTCRF01_AE': 0.4, '82VTCRF02_AG': 0.8, '82VTOther': 1.9, '82VTAll': 2.1, '82VFA': 1.8, '82VFB': 1.7, '82VFC': 0.6, '82VFD': 1.3, '82VFF': 1.4, '82VFG': 2.4, '82VFCRF01_AE': 3.0, '82VFCRF02_AG': 1.1, '82VFOther': 1.3, '82VFAll': 1.6, '82VSA': 0.0, '82VSB': 1.0, '82VSC': 0.1, '82VSD': 1.3, '82VSF': 0.5, '82VSG': 2.1, '82VSCRF01_AE': 0.1, '82VSCRF02_AG': 0.8, '82VSOther': 1.6, '82VSAll': 0.9, '82VMA': 0.0, '82VMB': 0.3, '82VMC': 0.5, '82VMD': 0.0, '82VMF': 0.5, '82VMG': 2.9, '82VMCRF01_AE': 0.0, '82VMCRF02_AG': 0.8, '82VMOther': 0.5, '82VMAll': 0.4, '82VCA': 0.0, '82VCB': 0.7, '82VCC': 0.4, '82VCD': 0.0, '82VCF': 0.1, '82VCG': 0.0, '82VCCRF01_AE': 0.0, '82VCCRF02_AG': 0.4, '82VCOther': 0.3, '82VCAll': 0.5, '82VLA': 0.0, '82VLB': 0.3, '82VLC': 0.3, '82VLD': 0.4, '82VLF': 0.4, '82VLG': 0.0, '82VLCRF01_AE': 0.0, '82VLCRF02_AG': 0.0, '82VLOther': 0.4, '82VLAll': 0.3, '82VNA': 0.0, '82VNB': 0.0, '82VNC': 0.0, '82VND': 0.0, '82VNF': 0.1, '82VNG': 0.3, '82VNCRF01_AE': 0.0, '82VNCRF02_AG': 0.0, '82VNOther': 0.0, '82VNAll': 0.0, '82VDA': 0.0, '82VDB': 0.0, '82VDC': 0.0, '82VDD': 0.0, '82VDF': 0.0, '82VDG': 0.0, '82VDCRF01_AE': 0.0, '82VDCRF02_AG': 0.4, '82VDOther': 0.0, '82VDAll': 0.0, '82VGA': 0.0, '82VGB': 0.0, '82VGC': 0.0, '82VGD': 0.0, '82VGF': 0.1, '82VGG': 0.0, '82VGCRF01_AE': 0.0, '82VGCRF02_AG': 0.0, '82VGOther': 0.2, '82VGAll': 0.0, '82VEA': 0.0, '82VEB': 0.0, '82VEC': 0.0, '82VED': 0.0, '82VEF': 0.0, '82VEG': 0.0, '82VECRF01_AE': 0.0, '82VECRF02_AG': 0.0, '82VEOther': 0.0, '82VEAll': 0.0, '82V~A': 0.0, '82V~B': 0.0, '82V~C': 0.0, '82V~D': 0.0, '82V~F': 0.0, '82V~G': 0.0, '82V~CRF01_AE': 0.0, '82V~CRF02_AG': 0.0, '82V~Other': 0.0, '82V~All': 0.0, '82VPA': 0.0, '82VPB': 0.0, '82VPC': 0.0, '82VPD': 0.0, '82VPF': 0.0, '82VPG': 0.0, '82VPCRF01_AE': 0.0, '82VPCRF02_AG': 0.0, '82VPOther': 0.0, '82VPAll': 0.0, '82VHA': 0.0, '82VHB': 0.0, '82VHC': 0.0, '82VHD': 0.0, '82VHF': 0.0, '82VHG': 0.0, '82VHCRF01_AE': 0.0, '82VHCRF02_AG': 0.0, '82VHOther': 0.0, '82VHAll': 0.0, '82VRA': 0.0, '82VRB': 0.0, '82VRC': 0.0, '82VRD': 0.0, '82VRF': 0.0, '82VRG': 0.0, '82VRCRF01_AE': 0.0, '82VRCRF02_AG': 0.0, '82VROther': 0.0, '82VRAll': 0.0, '83NDA': 3.6, '83NDB': 0.7, '83NDC': 0.4, '83NDD': 0.8, '83NDF': 2.3, '83NDG': 0.0, '83NDCRF01_AE': 4.5, '83NDCRF02_AG': 1.1, '83NDOther': 0.9, '83NDAll': 0.9, '83NSA': 0.0, '83NSB': 0.4, '83NSC': 0.1, '83NSD': 0.4, '83NSF': 0.0, '83NSG': 0.0, '83NSCRF01_AE': 0.1, '83NSCRF02_AG': 0.4, '83NSOther': 0.2, '83NSAll': 0.3, '83NHA': 0.0, '83NHB': 0.0, '83NHC': 0.0, '83NHD': 0.0, '83NHF': 0.0, '83NHG': 0.0, '83NHCRF01_AE': 0.0, '83NHCRF02_AG': 0.0, '83NHOther': 0.0, '83NHAll': 0.0, '83NYA': 0.0, '83NYB': 0.0, '83NYC': 0.0, '83NYD': 0.0, '83NYF': 0.0, '83NYG': 0.0, '83NYCRF01_AE': 0.0, '83NYCRF02_AG': 0.0, '83NYOther': 0.0, '83NYAll': 0.0, '83NPA': 0.0, '83NPB': 0.0, '83NPC': 0.0, '83NPD': 0.0, '83NPF': 0.0, '83NPG': 0.0, '83NPCRF01_AE': 0.0, '83NPCRF02_AG': 0.0, '83NPOther': 0.0, '83NPAll': 0.0, '83NTA': 0.0, '83NTB': 0.0, '83NTC': 0.0, '83NTD': 0.0, '83NTF': 0.0, '83NTG': 0.0, '83NTCRF01_AE': 0.0, '83NTCRF02_AG': 0.0, '83NTOther': 0.0, '83NTAll': 0.0, '83NKA': 0.0, '83NKB': 0.0, '83NKC': 0.0, '83NKD': 0.0, '83NKF': 0.0, '83NKG': 0.0, '83NKCRF01_AE': 0.0, '83NKCRF02_AG': 0.0, '83NKOther': 0.0, '83NKAll': 0.0, '83NLA': 0.0, '83NLB': 0.0, '83NLC': 0.0, '83NLD': 0.0, '83NLF': 0.0, '83NLG': 0.0, '83NLCRF01_AE': 0.0, '83NLCRF02_AG': 0.0, '83NLOther': 0.0, '83NLAll': 0.0, '83NIA': 0.0, '83NIB': 0.0, '83NIC': 0.0, '83NID': 0.0, '83NIF': 0.0, '83NIG': 0.0, '83NICRF01_AE': 0.0, '83NICRF02_AG': 0.0, '83NIOther': 0.0, '83NIAll': 0.0, '84IVA': 4.2, '84IVB': 13.5, '84IVC': 2.2, '84IVD': 6.3, '84IVF': 3.8, '84IVG': 6.2, '84IVCRF01_AE': 3.9, '84IVCRF02_AG': 13.8, '84IVOther': 3.5, '84IVAll': 10.2, '84ICA': 0.0, '84ICB': 0.1, '84ICC': 0.0, '84ICD': 0.4, '84ICF': 0.0, '84ICG': 0.0, '84ICCRF01_AE': 0.0, '84ICCRF02_AG': 0.4, '84ICOther': 0.1, '84ICAll': 0.1, '84IAA': 0.0, '84IAB': 0.1, '84IAC': 0.0, '84IAD': 0.0, '84IAF': 0.0, '84IAG': 0.0, '84IACRF01_AE': 0.0, '84IACRF02_AG': 0.4, '84IAOther': 0.0, '84IAAll': 0.1, '84IMA': 0.0, '84IMB': 0.0, '84IMC': 0.0, '84IMD': 0.0, '84IMF': 0.0, '84IMG': 0.0, '84IMCRF01_AE': 0.0, '84IMCRF02_AG': 0.0, '84IMOther': 0.0, '84IMAll': 0.0, '84ILA': 0.0, '84ILB': 0.0, '84ILC': 0.0, '84ILD': 0.0, '84ILF': 0.0, '84ILG': 0.0, '84ILCRF01_AE': 0.1, '84ILCRF02_AG': 0.0, '84ILOther': 0.0, '84ILAll': 0.0, '84ISA': 0.0, '84ISB': 0.0, '84ISC': 0.0, '84ISD': 0.0, '84ISF': 0.0, '84ISG': 0.0, '84ISCRF01_AE': 0.0, '84ISCRF02_AG': 0.0, '84ISOther': 0.0, '84ISAll': 0.0, '84ITA': 0.0, '84ITB': 0.0, '84ITC': 0.0, '84ITD': 0.0, '84ITF': 0.0, '84ITG': 0.0, '84ITCRF01_AE': 0.0, '84ITCRF02_AG': 0.0, '84ITOther': 0.0, '84ITAll': 0.0, '84IKA': 0.0, '84IKB': 0.0, '84IKC': 0.0, '84IKD': 0.0, '84IKF': 0.0, '84IKG': 0.0, '84IKCRF01_AE': 0.0, '84IKCRF02_AG': 0.0, '84IKOther': 0.0, '84IKAll': 0.0, '84IRA': 0.0, '84IRB': 0.0, '84IRC': 0.0, '84IRD': 0.0, '84IRF': 0.0, '84IRG': 0.0, '84IRCRF01_AE': 0.0, '84IRCRF02_AG': 0.0, '84IROther': 0.0, '84IRAll': 0.0, '85IVA': 0.9, '85IVB': 4.7, '85IVC': 2.3, '85IVD': 3.3, '85IVF': 3.9, '85IVG': 0.8, '85IVCRF01_AE': 0.3, '85IVCRF02_AG': 0.0, '85IVOther': 4.3, '85IVAll': 4.1, '85IMA': 0.0, '85IMB': 0.0, '85IMC': 0.0, '85IMD': 0.4, '85IMF': 0.0, '85IMG': 0.3, '85IMCRF01_AE': 0.0, '85IMCRF02_AG': 0.0, '85IMOther': 0.0, '85IMAll': 0.0, '85ILA': 0.0, '85ILB': 0.1, '85ILC': 0.0, '85ILD': 0.0, '85ILF': 0.0, '85ILG': 0.0, '85ILCRF01_AE': 0.0, '85ILCRF02_AG': 0.0, '85ILOther': 0.1, '85ILAll': 0.0, '85ITA': 0.0, '85ITB': 0.0, '85ITC': 0.0, '85ITD': 0.0, '85ITF': 0.0, '85ITG': 0.0, '85ITCRF01_AE': 0.0, '85ITCRF02_AG': 0.0, '85ITOther': 0.0, '85ITAll': 0.0, '85INA': 0.0, '85INB': 0.0, '85INC': 0.0, '85IND': 0.0, '85INF': 0.0, '85ING': 0.0, '85INCRF01_AE': 0.0, '85INCRF02_AG': 0.0, '85INOther': 0.0, '85INAll': 0.0, '85IFA': 0.0, '85IFB': 0.0, '85IFC': 0.0, '85IFD': 0.0, '85IFF': 0.0, '85IFG': 0.0, '85IFCRF01_AE': 0.0, '85IFCRF02_AG': 0.0, '85IFOther': 0.0, '85IFAll': 0.0, '85ISA': 0.0, '85ISB': 0.0, '85ISC': 0.0, '85ISD': 0.0, '85ISF': 0.0, '85ISG': 0.0, '85ISCRF01_AE': 0.0, '85ISCRF02_AG': 0.0, '85ISOther': 0.0, '85ISAll': 0.0, '85IKA': 0.0, '85IKB': 0.0, '85IKC': 0.0, '85IKD': 0.0, '85IKF': 0.0, '85IKG': 0.0, '85IKCRF01_AE': 0.0, '85IKCRF02_AG': 0.0, '85IKOther': 0.0, '85IKAll': 0.0, '85IRA': 0.0, '85IRB': 0.0, '85IRC': 0.0, '85IRD': 0.0, '85IRF': 0.0, '85IRG': 0.0, '85IRCRF01_AE': 0.0, '85IRCRF02_AG': 0.0, '85IROther': 0.0, '85IRAll': 0.0, '86GRA': 0.0, '86GRB': 0.0, '86GRC': 0.0, '86GRD': 0.0, '86GRF': 0.1, '86GRG': 0.0, '86GRCRF01_AE': 0.0, '86GRCRF02_AG': 0.0, '86GROther': 0.1, '86GRAll': 0.0, '86GEA': 0.0, '86GEB': 0.0, '86GEC': 0.0, '86GED': 0.0, '86GEF': 0.0, '86GEG': 0.0, '86GECRF01_AE': 0.0, '86GECRF02_AG': 0.0, '86GEOther': 0.1, '86GEAll': 0.0, '86GAA': 0.0, '86GAB': 0.0, '86GAC': 0.0, '86GAD': 0.0, '86GAF': 0.1, '86GAG': 0.0, '86GACRF01_AE': 0.0, '86GACRF02_AG': 0.0, '86GAOther': 0.0, '86GAAll': 0.0, '86GWA': 0.0, '86GWB': 0.0, '86GWC': 0.0, '86GWD': 0.0, '86GWF': 0.0, '86GWG': 0.0, '86GWCRF01_AE': 0.0, '86GWCRF02_AG': 0.0, '86GWOther': 0.0, '86GWAll': 0.0, '86GKA': 0.0, '86GKB': 0.0, '86GKC': 0.0, '86GKD': 0.0, '86GKF': 0.0, '86GKG': 0.0, '86GKCRF01_AE': 0.0, '86GKCRF02_AG': 0.0, '86GKOther': 0.0, '86GKAll': 0.0, '86GVA': 0.0, '86GVB': 0.0, '86GVC': 0.0, '86GVD': 0.0, '86GVF': 0.0, '86GVG': 0.0, '86GVCRF01_AE': 0.0, '86GVCRF02_AG': 0.0, '86GVOther': 0.0, '86GVAll': 0.0, '87RKA': 0.6, '87RKB': 0.1, '87RKC': 0.0, '87RKD': 0.0, '87RKF': 0.1, '87RKG': 0.0, '87RKCRF01_AE': 0.0, '87RKCRF02_AG': 0.0, '87RKOther': 0.0, '87RKAll': 0.1, '87RQA': 0.0, '87RQB': 0.0, '87RQC': 0.0, '87RQD': 0.0, '87RQF': 0.0, '87RQG': 0.0, '87RQCRF01_AE': 0.6, '87RQCRF02_AG': 0.0, '87RQOther': 0.0, '87RQAll': 0.0, '87RGA': 0.0, '87RGB': 0.0, '87RGC': 0.0, '87RGD': 0.0, '87RGF': 0.0, '87RGG': 0.0, '87RGCRF01_AE': 0.0, '87RGCRF02_AG': 0.4, '87RGOther': 0.1, '87RGAll': 0.0, '87REA': 0.0, '87REB': 0.0, '87REC': 0.0, '87RED': 0.0, '87REF': 0.0, '87REG': 0.3, '87RECRF01_AE': 0.0, '87RECRF02_AG': 0.0, '87REOther': 0.0, '87REAll': 0.0, '87RTA': 0.0, '87RTB': 0.1, '87RTC': 0.0, '87RTD': 0.0, '87RTF': 0.0, '87RTG': 0.0, '87RTCRF01_AE': 0.0, '87RTCRF02_AG': 0.0, '87RTOther': 0.0, '87RTAll': 0.1, '87RIA': 0.0, '87RIB': 0.1, '87RIC': 0.0, '87RID': 0.0, '87RIF': 0.1, '87RIG': 0.0, '87RICRF01_AE': 0.0, '87RICRF02_AG': 0.0, '87RIOther': 0.0, '87RIAll': 0.0, '87RSA': 0.0, '87RSB': 0.0, '87RSC': 0.0, '87RSD': 0.0, '87RSF': 0.0, '87RSG': 0.0, '87RSCRF01_AE': 0.0, '87RSCRF02_AG': 0.0, '87RSOther': 0.0, '87RSAll': 0.0, '87RVA': 0.0, '87RVB': 0.0, '87RVC': 0.0, '87RVD': 0.0, '87RVF': 0.0, '87RVG': 0.0, '87RVCRF01_AE': 0.0, '87RVCRF02_AG': 0.0, '87RVOther': 0.0, '87RVAll': 0.0, '87RLA': 0.0, '87RLB': 0.0, '87RLC': 0.0, '87RLD': 0.0, '87RLF': 0.0, '87RLG': 0.0, '87RLCRF01_AE': 0.0, '87RLCRF02_AG': 0.0, '87RLOther': 0.0, '87RLAll': 0.0, '87RPA': 0.0, '87RPB': 0.0, '87RPC': 0.0, '87RPD': 0.0, '87RPF': 0.0, '87RPG': 0.0, '87RPCRF01_AE': 0.0, '87RPCRF02_AG': 0.0, '87RPOther': 0.0, '87RPAll': 0.0, '87RHA': 0.0, '87RHB': 0.0, '87RHC': 0.0, '87RHD': 0.0, '87RHF': 0.0, '87RHG': 0.0, '87RHCRF01_AE': 0.0, '87RHCRF02_AG': 0.0, '87RHOther': 0.0, '87RHAll': 0.0, '88NDA': 0.0, '88NDB': 5.7, '88NDC': 1.8, '88NDD': 4.6, '88NDF': 4.1, '88NDG': 1.0, '88NDCRF01_AE': 0.6, '88NDCRF02_AG': 0.0, '88NDOther': 5.7, '88NDAll': 4.8, '88NSA': 2.7, '88NSB': 1.4, '88NSC': 0.9, '88NSD': 0.8, '88NSF': 2.8, '88NSG': 1.6, '88NSCRF01_AE': 1.9, '88NSCRF02_AG': 2.3, '88NSOther': 1.3, '88NSAll': 1.5, '88NGA': 0.9, '88NGB': 0.1, '88NGC': 0.1, '88NGD': 0.8, '88NGF': 0.3, '88NGG': 0.3, '88NGCRF01_AE': 0.1, '88NGCRF02_AG': 0.0, '88NGOther': 0.1, '88NGAll': 0.1, '88NKA': 0.0, '88NKB': 0.0, '88NKC': 0.0, '88NKD': 0.4, '88NKF': 0.0, '88NKG': 0.0, '88NKCRF01_AE': 0.0, '88NKCRF02_AG': 0.4, '88NKOther': 0.0, '88NKAll': 0.0, '88NTA': 0.0, '88NTB': 0.1, '88NTC': 0.1, '88NTD': 0.0, '88NTF': 0.3, '88NTG': 0.0, '88NTCRF01_AE': 0.0, '88NTCRF02_AG': 0.0, '88NTOther': 0.1, '88NTAll': 0.1, '88NHA': 0.0, '88NHB': 0.0, '88NHC': 0.0, '88NHD': 0.0, '88NHF': 0.0, '88NHG': 0.0, '88NHCRF01_AE': 0.1, '88NHCRF02_AG': 0.0, '88NHOther': 0.0, '88NHAll': 0.0, '88NEA': 0.0, '88NEB': 0.0, '88NEC': 0.0, '88NED': 0.0, '88NEF': 0.1, '88NEG': 0.0, '88NECRF01_AE': 0.0, '88NECRF02_AG': 0.0, '88NEOther': 0.0, '88NEAll': 0.0, '88NYA': 0.0, '88NYB': 0.0, '88NYC': 0.0, '88NYD': 0.0, '88NYF': 0.0, '88NYG': 0.0, '88NYCRF01_AE': 0.0, '88NYCRF02_AG': 0.0, '88NYOther': 0.1, '88NYAll': 0.0, '88NIA': 0.0, '88NIB': 0.0, '88NIC': 0.0, '88NID': 0.0, '88NIF': 0.0, '88NIG': 0.0, '88NICRF01_AE': 0.0, '88NICRF02_AG': 0.0, '88NIOther': 0.0, '88NIAll': 0.0, '89LMA': 80.6, '89LMB': 2.7, '89LMC': 72.9, '89LMD': 3.4, '89LMF': 53.5, '89LMG': 63.3, '89LMCRF01_AE': 84.0, '89LMCRF02_AG': 71.9, '89LMOther': 12.7, '89LMAll': 19.7, '89LIA': 11.6, '89LIB': 1.0, '89LIC': 8.9, '89LID': 1.3, '89LIF': 10.0, '89LIG': 22.3, '89LICRF01_AE': 8.1, '89LICRF02_AG': 13.8, '89LIOther': 2.3, '89LIAll': 3.3, '89LVA': 4.8, '89LVB': 3.6, '89LVC': 1.1, '89LVD': 5.9, '89LVF': 1.7, '89LVG': 12.1, '89LVCRF01_AE': 4.0, '89LVCRF02_AG': 10.8, '89LVOther': 1.5, '89LVAll': 3.2, '89LTA': 2.1, '89LTB': 0.2, '89LTC': 0.8, '89LTD': 0.0, '89LTF': 2.5, '89LTG': 0.5, '89LTCRF01_AE': 1.0, '89LTCRF02_AG': 2.3, '89LTOther': 0.4, '89LTAll': 0.5, '89LFA': 0.0, '89LFB': 0.3, '89LFC': 0.2, '89LFD': 1.3, '89LFF': 0.3, '89LFG': 0.0, '89LFCRF01_AE': 0.0, '89LFCRF02_AG': 0.0, '89LFOther': 0.1, '89LFAll': 0.2, '89LPA': 0.0, '89LPB': 0.0, '89LPC': 0.0, '89LPD': 0.0, '89LPF': 0.1, '89LPG': 0.0, '89LPCRF01_AE': 0.0, '89LPCRF02_AG': 0.0, '89LPOther': 0.2, '89LPAll': 0.0, '89LKA': 0.0, '89LKB': 0.0, '89LKC': 0.0, '89LKD': 0.0, '89LKF': 0.0, '89LKG': 0.3, '89LKCRF01_AE': 0.0, '89LKCRF02_AG': 0.0, '89LKOther': 0.0, '89LKAll': 0.0, '89LRA': 0.0, '89LRB': 0.0, '89LRC': 0.0, '89LRD': 0.0, '89LRF': 0.0, '89LRG': 0.0, '89LRCRF01_AE': 0.0, '89LRCRF02_AG': 0.0, '89LROther': 0.1, '89LRAll': 0.0, '89LSA': 0.0, '89LSB': 0.0, '89LSC': 0.0, '89LSD': 0.0, '89LSF': 0.0, '89LSG': 0.0, '89LSCRF01_AE': 0.0, '89LSCRF02_AG': 0.0, '89LSOther': 0.0, '89LSAll': 0.0, '89LYA': 0.0, '89LYB': 0.0, '89LYC': 0.0, '89LYD': 0.0, '89LYF': 0.0, '89LYG': 0.0, '89LYCRF01_AE': 0.0, '89LYCRF02_AG': 0.0, '89LYOther': 0.0, '89LYAll': 0.0, '89LQA': 0.0, '89LQB': 0.0, '89LQC': 0.0, '89LQD': 0.0, '89LQF': 0.0, '89LQG': 0.0, '89LQCRF01_AE': 0.0, '89LQCRF02_AG': 0.0, '89LQOther': 0.0, '89LQAll': 0.0, '89LAA': 0.0, '89LAB': 0.0, '89LAC': 0.0, '89LAD': 0.0, '89LAF': 0.0, '89LAG': 0.0, '89LACRF01_AE': 0.0, '89LACRF02_AG': 0.0, '89LAOther': 0.0, '89LAAll': 0.0, '90LMA': 7.2, '90LMB': 30.6, '90LMC': 8.4, '90LMD': 20.5, '90LMF': 18.3, '90LMG': 31.5, '90LMCRF01_AE': 6.1, '90LMCRF02_AG': 8.5, '90LMOther': 21.4, '90LMAll': 25.3, '90LVA': 0.0, '90LVB': 0.0, '90LVC': 0.0, '90LVD': 0.4, '90LVF': 0.0, '90LVG': 0.0, '90LVCRF01_AE': 0.1, '90LVCRF02_AG': 0.0, '90LVOther': 0.0, '90LVAll': 0.0, '90LFA': 0.0, '90LFB': 0.0, '90LFC': 0.0, '90LFD': 0.0, '90LFF': 0.0, '90LFG': 0.0, '90LFCRF01_AE': 0.0, '90LFCRF02_AG': 0.0, '90LFOther': 0.1, '90LFAll': 0.0, '90LSA': 0.0, '90LSB': 0.0, '90LSC': 0.0, '90LSD': 0.0, '90LSF': 0.0, '90LSG': 0.0, '90LSCRF01_AE': 0.0, '90LSCRF02_AG': 0.0, '90LSOther': 0.0, '90LSAll': 0.0, '90LWA': 0.0, '90LWB': 0.0, '90LWC': 0.0, '90LWD': 0.0, '90LWF': 0.0, '90LWG': 0.0, '90LWCRF01_AE': 0.1, '90LWCRF02_AG': 0.0, '90LWOther': 0.0, '90LWAll': 0.0, '90LRA': 0.0, '90LRB': 0.0, '90LRC': 0.0, '90LRD': 0.0, '90LRF': 0.0, '90LRG': 0.0, '90LRCRF01_AE': 0.0, '90LRCRF02_AG': 0.0, '90LROther': 0.1, '90LRAll': 0.0, '90LKA': 0.0, '90LKB': 0.0, '90LKC': 0.0, '90LKD': 0.0, '90LKF': 0.0, '90LKG': 0.0, '90LKCRF01_AE': 0.0, '90LKCRF02_AG': 0.0, '90LKOther': 0.0, '90LKAll': 0.0, '90LCA': 0.0, '90LCB': 0.0, '90LCC': 0.0, '90LCD': 0.0, '90LCF': 0.0, '90LCG': 0.0, '90LCCRF01_AE': 0.0, '90LCCRF02_AG': 0.0, '90LCOther': 0.0, '90LCAll': 0.0, '90LPA': 0.0, '90LPB': 0.0, '90LPC': 0.0, '90LPD': 0.0, '90LPF': 0.0, '90LPG': 0.0, '90LPCRF01_AE': 0.0, '90LPCRF02_AG': 0.0, '90LPOther': 0.0, '90LPAll': 0.0, '90LIA': 0.0, '90LIB': 0.0, '90LIC': 0.0, '90LID': 0.0, '90LIF': 0.0, '90LIG': 0.0, '90LICRF01_AE': 0.0, '90LICRF02_AG': 0.0, '90LIOther': 0.0, '90LIAll': 0.0, '91TSA': 1.2, '91TSB': 1.9, '91TSC': 0.4, '91TSD': 1.7, '91TSF': 1.2, '91TSG': 2.1, '91TSCRF01_AE': 0.3, '91TSCRF02_AG': 3.4, '91TSOther': 0.8, '91TSAll': 1.5, '91TVA': 4.5, '91TVB': 0.0, '91TVC': 0.0, '91TVD': 0.0, '91TVF': 0.1, '91TVG': 0.0, '91TVCRF01_AE': 0.3, '91TVCRF02_AG': 0.0, '91TVOther': 0.1, '91TVAll': 0.1, '91TAA': 0.6, '91TAB': 0.5, '91TAC': 0.1, '91TAD': 0.0, '91TAF': 0.2, '91TAG': 1.3, '91TACRF01_AE': 0.4, '91TACRF02_AG': 1.1, '91TAOther': 0.3, '91TAAll': 0.4, '91TIA': 0.0, '91TIB': 0.0, '91TIC': 0.2, '91TID': 0.0, '91TIF': 0.1, '91TIG': 0.0, '91TICRF01_AE': 0.3, '91TICRF02_AG': 0.0, '91TIOther': 0.1, '91TIAll': 0.1, '91TNA': 0.0, '91TNB': 0.1, '91TNC': 0.1, '91TND': 0.0, '91TNF': 0.0, '91TNG': 0.5, '91TNCRF01_AE': 0.0, '91TNCRF02_AG': 0.0, '91TNOther': 0.0, '91TNAll': 0.1, '91TPA': 0.0, '91TPB': 0.0, '91TPC': 0.0, '91TPD': 0.0, '91TPF': 0.0, '91TPG': 0.0, '91TPCRF01_AE': 0.1, '91TPCRF02_AG': 0.0, '91TPOther': 0.0, '91TPAll': 0.0, '91TCA': 0.0, '91TCB': 0.0, '91TCC': 0.0, '91TCD': 0.0, '91TCF': 0.0, '91TCG': 0.3, '91TCCRF01_AE': 0.0, '91TCCRF02_AG': 0.0, '91TCOther': 0.0, '91TCAll': 0.0, '91TGA': 0.0, '91TGB': 0.0, '91TGC': 0.0, '91TGD': 0.0, '91TGF': 0.0, '91TGG': 0.0, '91TGCRF01_AE': 0.1, '91TGCRF02_AG': 0.0, '91TGOther': 0.0, '91TGAll': 0.0, '91TKA': 0.0, '91TKB': 0.0, '91TKC': 0.0, '91TKD': 0.0, '91TKF': 0.0, '91TKG': 0.0, '91TKCRF01_AE': 0.0, '91TKCRF02_AG': 0.0, '91TKOther': 0.0, '91TKAll': 0.0, '91TMA': 0.0, '91TMB': 0.0, '91TMC': 0.0, '91TMD': 0.0, '91TMF': 0.0, '91TMG': 0.0, '91TMCRF01_AE': 0.0, '91TMCRF02_AG': 0.0, '91TMOther': 0.0, '91TMAll': 0.0, '91TDA': 0.0, '91TDB': 0.0, '91TDC': 0.0, '91TDD': 0.0, '91TDF': 0.0, '91TDG': 0.0, '91TDCRF01_AE': 0.0, '91TDCRF02_AG': 0.0, '91TDOther': 0.0, '91TDAll': 0.0, '92QKA': 0.6, '92QKB': 3.2, '92QKC': 0.6, '92QKD': 3.3, '92QKF': 5.0, '92QKG': 2.3, '92QKCRF01_AE': 0.7, '92QKCRF02_AG': 1.5, '92QKOther': 2.7, '92QKAll': 2.9, '92QRA': 0.0, '92QRB': 0.9, '92QRC': 0.1, '92QRD': 0.0, '92QRF': 1.0, '92QRG': 0.8, '92QRCRF01_AE': 0.0, '92QRCRF02_AG': 0.0, '92QROther': 0.6, '92QRAll': 0.7, '92QEA': 0.0, '92QEB': 0.3, '92QEC': 0.0, '92QED': 0.0, '92QEF': 0.0, '92QEG': 0.0, '92QECRF01_AE': 0.3, '92QECRF02_AG': 0.0, '92QEOther': 0.2, '92QEAll': 0.2, '92QHA': 0.0, '92QHB': 0.1, '92QHC': 0.0, '92QHD': 0.0, '92QHF': 0.1, '92QHG': 0.0, '92QHCRF01_AE': 0.0, '92QHCRF02_AG': 0.0, '92QHOther': 0.0, '92QHAll': 0.1, '92QLA': 0.3, '92QLB': 0.0, '92QLC': 0.0, '92QLD': 0.0, '92QLF': 0.0, '92QLG': 0.0, '92QLCRF01_AE': 0.0, '92QLCRF02_AG': 0.0, '92QLOther': 0.0, '92QLAll': 0.0, '92QTA': 0.0, '92QTB': 0.0, '92QTC': 0.0, '92QTD': 0.0, '92QTF': 0.1, '92QTG': 0.0, '92QTCRF01_AE': 0.0, '92QTCRF02_AG': 0.0, '92QTOther': 0.0, '92QTAll': 0.0, '92QMA': 0.0, '92QMB': 0.0, '92QMC': 0.0, '92QMD': 0.0, '92QMF': 0.1, '92QMG': 0.0, '92QMCRF01_AE': 0.0, '92QMCRF02_AG': 0.0, '92QMOther': 0.0, '92QMAll': 0.0, '92QPA': 0.0, '92QPB': 0.0, '92QPC': 0.0, '92QPD': 0.0, '92QPF': 0.0, '92QPG': 0.0, '92QPCRF01_AE': 0.0, '92QPCRF02_AG': 0.0, '92QPOther': 0.0, '92QPAll': 0.0, '92QSA': 0.0, '92QSB': 0.0, '92QSC': 0.0, '92QSD': 0.0, '92QSF': 0.0, '92QSG': 0.0, '92QSCRF01_AE': 0.0, '92QSCRF02_AG': 0.0, '92QSOther': 0.0, '92QSAll': 0.0, '92QNA': 0.0, '92QNB': 0.0, '92QNC': 0.0, '92QND': 0.0, '92QNF': 0.0, '92QNG': 0.0, '92QNCRF01_AE': 0.0, '92QNCRF02_AG': 0.0, '92QNOther': 0.0, '92QNAll': 0.0, '92QAA': 0.0, '92QAB': 0.0, '92QAC': 0.0, '92QAD': 0.0, '92QAF': 0.0, '92QAG': 0.0, '92QACRF01_AE': 0.0, '92QACRF02_AG': 0.0, '92QAOther': 0.0, '92QAAll': 0.0, '92QWA': 0.0, '92QWB': 0.0, '92QWC': 0.0, '92QWD': 0.0, '92QWF': 0.0, '92QWG': 0.0, '92QWCRF01_AE': 0.0, '92QWCRF02_AG': 0.0, '92QWOther': 0.0, '92QWAll': 0.0, '92QDA': 0.0, '92QDB': 0.0, '92QDC': 0.0, '92QDD': 0.0, '92QDF': 0.0, '92QDG': 0.0, '92QDCRF01_AE': 0.0, '92QDCRF02_AG': 0.0, '92QDOther': 0.0, '92QDAll': 0.0, '93ILA': 21.2, '93ILB': 36.7, '93ILC': 97.1, '93ILD': 21.4, '93ILF': 18.1, '93ILG': 2.3, '93ILCRF01_AE': 12.0, '93ILCRF02_AG': 1.5, '93ILOther': 48.4, '93ILAll': 41.7, '93IMA': 1.5, '93IMB': 0.6, '93IMC': 0.0, '93IMD': 0.4, '93IMF': 0.2, '93IMG': 2.9, '93IMCRF01_AE': 1.0, '93IMCRF02_AG': 1.1, '93IMOther': 0.2, '93IMAll': 0.5, '93IVA': 0.9, '93IVB': 0.1, '93IVC': 0.0, '93IVD': 0.0, '93IVF': 0.1, '93IVG': 0.3, '93IVCRF01_AE': 0.9, '93IVCRF02_AG': 0.8, '93IVOther': 0.1, '93IVAll': 0.1, '93IFA': 0.0, '93IFB': 0.0, '93IFC': 0.0, '93IFD': 0.0, '93IFF': 0.1, '93IFG': 0.0, '93IFCRF01_AE': 0.1, '93IFCRF02_AG': 0.0, '93IFOther': 0.1, '93IFAll': 0.0, '93IDA': 0.0, '93IDB': 0.0, '93IDC': 0.0, '93IDD': 0.0, '93IDF': 0.0, '93IDG': 0.0, '93IDCRF01_AE': 0.0, '93IDCRF02_AG': 0.0, '93IDOther': 0.0, '93IDAll': 0.0, '93ISA': 0.0, '93ISB': 0.0, '93ISC': 0.0, '93ISD': 0.0, '93ISF': 0.1, '93ISG': 0.0, '93ISCRF01_AE': 0.0, '93ISCRF02_AG': 0.0, '93ISOther': 0.0, '93ISAll': 0.0, '93ITA': 0.0, '93ITB': 0.0, '93ITC': 0.0, '93ITD': 0.0, '93ITF': 0.0, '93ITG': 0.0, '93ITCRF01_AE': 0.0, '93ITCRF02_AG': 0.0, '93ITOther': 0.0, '93ITAll': 0.0, '93INA': 0.0, '93INB': 0.0, '93INC': 0.0, '93IND': 0.0, '93INF': 0.0, '93ING': 0.0, '93INCRF01_AE': 0.0, '93INCRF02_AG': 0.0, '93INOther': 0.0, '93INAll': 0.0, '93IYA': 0.0, '93IYB': 0.0, '93IYC': 0.0, '93IYD': 0.0, '93IYF': 0.0, '93IYG': 0.0, '93IYCRF01_AE': 0.0, '93IYCRF02_AG': 0.0, '93IYOther': 0.0, '93IYAll': 0.0, '93IAA': 0.0, '93IAB': 0.0, '93IAC': 0.0, '93IAD': 0.0, '93IAF': 0.0, '93IAG': 0.0, '93IACRF01_AE': 0.0, '93IACRF02_AG': 0.0, '93IAOther': 0.0, '93IAAll': 0.0, '93IPA': 0.0, '93IPB': 0.0, '93IPC': 0.0, '93IPD': 0.0, '93IPF': 0.0, '93IPG': 0.0, '93IPCRF01_AE': 0.0, '93IPCRF02_AG': 0.0, '93IPOther': 0.0, '93IPAll': 0.0, '93IHA': 0.0, '93IHB': 0.0, '93IHC': 0.0, '93IHD': 0.0, '93IHF': 0.0, '93IHG': 0.0, '93IHCRF01_AE': 0.0, '93IHCRF02_AG': 0.0, '93IHOther': 0.0, '93IHAll': 0.0, '94GAA': 0.3, '94GAB': 0.0, '94GAC': 0.0, '94GAD': 0.0, '94GAF': 0.0, '94GAG': 0.0, '94GACRF01_AE': 0.1, '94GACRF02_AG': 0.0, '94GAOther': 0.0, '94GAAll': 0.0, '94GVA': 0.0, '94GVB': 0.0, '94GVC': 0.0, '94GVD': 0.0, '94GVF': 0.0, '94GVG': 0.0, '94GVCRF01_AE': 0.0, '94GVCRF02_AG': 0.0, '94GVOther': 0.0, '94GVAll': 0.0, '94GNA': 0.0, '94GNB': 0.0, '94GNC': 0.0, '94GND': 0.0, '94GNF': 0.0, '94GNG': 0.0, '94GNCRF01_AE': 0.0, '94GNCRF02_AG': 0.0, '94GNOther': 0.0, '94GNAll': 0.0, '94GCA': 0.0, '94GCB': 0.0, '94GCC': 0.0, '94GCD': 0.0, '94GCF': 0.0, '94GCG': 0.0, '94GCCRF01_AE': 0.0, '94GCCRF02_AG': 0.0, '94GCOther': 0.0, '94GCAll': 0.0, '94GDA': 0.0, '94GDB': 0.0, '94GDC': 0.0, '94GDD': 0.0, '94GDF': 0.1, '94GDG': 0.0, '94GDCRF01_AE': 0.0, '94GDCRF02_AG': 0.0, '94GDOther': 0.0, '94GDAll': 0.0, '94GSA': 0.0, '94GSB': 0.0, '94GSC': 0.0, '94GSD': 0.0, '94GSF': 0.0, '94GSG': 0.0, '94GSCRF01_AE': 0.0, '94GSCRF02_AG': 0.0, '94GSOther': 0.0, '94GSAll': 0.0, '94GEA': 0.0, '94GEB': 0.0, '94GEC': 0.0, '94GED': 0.0, '94GEF': 0.0, '94GEG': 0.0, '94GECRF01_AE': 0.0, '94GECRF02_AG': 0.0, '94GEOther': 0.0, '94GEAll': 0.0, '94G~A': 0.0, '94G~B': 0.0, '94G~C': 0.0, '94G~D': 0.0, '94G~F': 0.0, '94G~G': 0.0, '94G~CRF01_AE': 0.0, '94G~CRF02_AG': 0.0, '94G~Other': 0.0, '94G~All': 0.0, '94GQA': 0.0, '94GQB': 0.0, '94GQC': 0.0, '94GQD': 0.0, '94GQF': 0.0, '94GQG': 0.0, '94GQCRF01_AE': 0.0, '94GQCRF02_AG': 0.0, '94GQOther': 0.0, '94GQAll': 0.0, '94GWA': 0.0, '94GWB': 0.0, '94GWC': 0.0, '94GWD': 0.0, '94GWF': 0.0, '94GWG': 0.0, '94GWCRF01_AE': 0.0, '94GWCRF02_AG': 0.0, '94GWOther': 0.0, '94GWAll': 0.0, '94GHA': 0.0, '94GHB': 0.0, '94GHC': 0.0, '94GHD': 0.0, '94GHF': 0.0, '94GHG': 0.0, '94GHCRF01_AE': 0.0, '94GHCRF02_AG': 0.0, '94GHOther': 0.0, '94GHAll': 0.0, '94GRA': 0.0, '94GRB': 0.0, '94GRC': 0.0, '94GRD': 0.0, '94GRF': 0.0, '94GRG': 0.0, '94GRCRF01_AE': 0.0, '94GRCRF02_AG': 0.0, '94GROther': 0.0, '94GRAll': 0.0, '95CFA': 0.0, '95CFB': 1.8, '95CFC': 0.3, '95CFD': 0.0, '95CFF': 0.3, '95CFG': 0.3, '95CFCRF01_AE': 0.0, '95CFCRF02_AG': 0.0, '95CFOther': 0.3, '95CFAll': 1.3, '95CSA': 0.0, '95CSB': 0.1, '95CSC': 0.0, '95CSD': 0.0, '95CSF': 0.1, '95CSG': 0.0, '95CSCRF01_AE': 0.0, '95CSCRF02_AG': 0.0, '95CSOther': 0.1, '95CSAll': 0.1, '95CVA': 0.0, '95CVB': 0.1, '95CVC': 0.0, '95CVD': 0.4, '95CVF': 0.1, '95CVG': 0.0, '95CVCRF01_AE': 0.0, '95CVCRF02_AG': 0.0, '95CVOther': 0.1, '95CVAll': 0.1, '95CLA': 0.0, '95CLB': 0.1, '95CLC': 0.0, '95CLD': 0.0, '95CLF': 0.2, '95CLG': 0.0, '95CLCRF01_AE': 0.0, '95CLCRF02_AG': 0.0, '95CLOther': 0.1, '95CLAll': 0.1, '95CWA': 0.0, '95CWB': 0.0, '95CWC': 0.0, '95CWD': 0.0, '95CWF': 0.0, '95CWG': 0.0, '95CWCRF01_AE': 0.0, '95CWCRF02_AG': 0.0, '95CWOther': 0.0, '95CWAll': 0.0, '95CGA': 0.0, '95CGB': 0.0, '95CGC': 0.0, '95CGD': 0.0, '95CGF': 0.0, '95CGG': 0.0, '95CGCRF01_AE': 0.0, '95CGCRF02_AG': 0.0, '95CGOther': 0.1, '95CGAll': 0.0, '95CNA': 0.0, '95CNB': 0.0, '95CNC': 0.0, '95CND': 0.0, '95CNF': 0.1, '95CNG': 0.0, '95CNCRF01_AE': 0.0, '95CNCRF02_AG': 0.0, '95CNOther': 0.0, '95CNAll': 0.0, '95CYA': 0.0, '95CYB': 0.0, '95CYC': 0.1, '95CYD': 0.0, '95CYF': 0.0, '95CYG': 0.0, '95CYCRF01_AE': 0.0, '95CYCRF02_AG': 0.0, '95CYOther': 0.0, '95CYAll': 0.0, '95CMA': 0.0, '95CMB': 0.0, '95CMC': 0.0, '95CMD': 0.0, '95CMF': 0.0, '95CMG': 0.0, '95CMCRF01_AE': 0.0, '95CMCRF02_AG': 0.0, '95CMOther': 0.0, '95CMAll': 0.0, '95CAA': 0.0, '95CAB': 0.0, '95CAC': 0.0, '95CAD': 0.0, '95CAF': 0.0, '95CAG': 0.0, '95CACRF01_AE': 0.0, '95CACRF02_AG': 0.0, '95CAOther': 0.0, '95CAAll': 0.0, '95CDA': 0.0, '95CDB': 0.0, '95CDC': 0.0, '95CDD': 0.0, '95CDF': 0.0, '95CDG': 0.0, '95CDCRF01_AE': 0.0, '95CDCRF02_AG': 0.0, '95CDOther': 0.0, '95CDAll': 0.0, '95CIA': 0.0, '95CIB': 0.0, '95CIC': 0.0, '95CID': 0.0, '95CIF': 0.0, '95CIG': 0.0, '95CICRF01_AE': 0.0, '95CICRF02_AG': 0.0, '95CIOther': 0.0, '95CIAll': 0.0, '95CRA': 0.0, '95CRB': 0.0, '95CRC': 0.0, '95CRD': 0.0, '95CRF': 0.0, '95CRG': 0.0, '95CRCRF01_AE': 0.0, '95CRCRF02_AG': 0.0, '95CROther': 0.0, '95CRAll': 0.0, '96TSA': 0.0, '96TSB': 0.2, '96TSC': 0.1, '96TSD': 0.0, '96TSF': 0.1, '96TSG': 0.0, '96TSCRF01_AE': 0.0, '96TSCRF02_AG': 0.0, '96TSOther': 0.6, '96TSAll': 0.2, '96TNA': 0.0, '96TNB': 0.0, '96TNC': 0.0, '96TND': 0.0, '96TNF': 0.1, '96TNG': 0.5, '96TNCRF01_AE': 0.2, '96TNCRF02_AG': 0.0, '96TNOther': 0.0, '96TNAll': 0.0, '96TPA': 0.3, '96TPB': 0.0, '96TPC': 0.0, '96TPD': 0.0, '96TPF': 0.0, '96TPG': 0.0, '96TPCRF01_AE': 0.2, '96TPCRF02_AG': 0.0, '96TPOther': 0.1, '96TPAll': 0.0, '96TIA': 0.0, '96TIB': 0.0, '96TIC': 0.0, '96TID': 0.0, '96TIF': 0.0, '96TIG': 0.3, '96TICRF01_AE': 0.0, '96TICRF02_AG': 0.0, '96TIOther': 0.1, '96TIAll': 0.0, '96TAA': 0.0, '96TAB': 0.0, '96TAC': 0.0, '96TAD': 0.0, '96TAF': 0.0, '96TAG': 0.3, '96TACRF01_AE': 0.0, '96TACRF02_AG': 0.0, '96TAOther': 0.1, '96TAAll': 0.0, '96THA': 0.3, '96THB': 0.0, '96THC': 0.0, '96THD': 0.0, '96THF': 0.0, '96THG': 0.0, '96THCRF01_AE': 0.0, '96THCRF02_AG': 0.0, '96THOther': 0.0, '96THAll': 0.0, '96TYA': 0.0, '96TYB': 0.0, '96TYC': 0.0, '96TYD': 0.0, '96TYF': 0.0, '96TYG': 0.0, '96TYCRF01_AE': 0.0, '96TYCRF02_AG': 0.0, '96TYOther': 0.0, '96TYAll': 0.0, '96TLA': 0.0, '96TLB': 0.0, '96TLC': 0.0, '96TLD': 0.0, '96TLF': 0.0, '96TLG': 0.0, '96TLCRF01_AE': 0.0, '96TLCRF02_AG': 0.0, '96TLOther': 0.0, '96TLAll': 0.0, '96TDA': 0.0, '96TDB': 0.0, '96TDC': 0.0, '96TDD': 0.0, '96TDF': 0.0, '96TDG': 0.0, '96TDCRF01_AE': 0.0, '96TDCRF02_AG': 0.0, '96TDOther': 0.0, '96TDAll': 0.0, '97LIA': 0.0, '97LIB': 0.1, '97LIC': 0.2, '97LID': 0.4, '97LIF': 0.1, '97LIG': 0.0, '97LICRF01_AE': 0.6, '97LICRF02_AG': 0.0, '97LIOther': 0.2, '97LIAll': 0.1, '97LFA': 0.0, '97LFB': 0.1, '97LFC': 0.0, '97LFD': 0.0, '97LFF': 0.0, '97LFG': 0.0, '97LFCRF01_AE': 0.0, '97LFCRF02_AG': 0.0, '97LFOther': 0.0, '97LFAll': 0.0, '97LMA': 0.0, '97LMB': 0.0, '97LMC': 0.0, '97LMD': 0.0, '97LMF': 0.0, '97LMG': 0.0, '97LMCRF01_AE': 0.2, '97LMCRF02_AG': 0.0, '97LMOther': 0.1, '97LMAll': 0.0, '97LSA': 0.3, '97LSB': 0.0, '97LSC': 0.0, '97LSD': 0.0, '97LSF': 0.0, '97LSG': 0.0, '97LSCRF01_AE': 0.0, '97LSCRF02_AG': 0.0, '97LSOther': 0.0, '97LSAll': 0.0, '97LTA': 0.0, '97LTB': 0.0, '97LTC': 0.0, '97LTD': 0.0, '97LTF': 0.0, '97LTG': 0.0, '97LTCRF01_AE': 0.0, '97LTCRF02_AG': 0.0, '97LTOther': 0.0, '97LTAll': 0.0, '97LKA': 0.0, '97LKB': 0.0, '97LKC': 0.0, '97LKD': 0.0, '97LKF': 0.0, '97LKG': 0.0, '97LKCRF01_AE': 0.0, '97LKCRF02_AG': 0.0, '97LKOther': 0.0, '97LKAll': 0.0, '97LYA': 0.0, '97LYB': 0.0, '97LYC': 0.0, '97LYD': 0.0, '97LYF': 0.0, '97LYG': 0.0, '97LYCRF01_AE': 0.0, '97LYCRF02_AG': 0.0, '97LYOther': 0.0, '97LYAll': 0.0, '97LVA': 0.0, '97LVB': 0.0, '97LVC': 0.0, '97LVD': 0.0, '97LVF': 0.0, '97LVG': 0.0, '97LVCRF01_AE': 0.0, '97LVCRF02_AG': 0.0, '97LVOther': 0.0, '97LVAll': 0.0, '97LQA': 0.0, '97LQB': 0.0, '97LQC': 0.0, '97LQD': 0.0, '97LQF': 0.0, '97LQG': 0.0, '97LQCRF01_AE': 0.0, '97LQCRF02_AG': 0.0, '97LQOther': 0.0, '97LQAll': 0.0, '97LPA': 0.0, '97LPB': 0.0, '97LPC': 0.0, '97LPD': 0.0, '97LPF': 0.0, '97LPG': 0.0, '97LPCRF01_AE': 0.0, '97LPCRF02_AG': 0.0, '97LPOther': 0.0, '97LPAll': 0.0, '97LRA': 0.0, '97LRB': 0.0, '97LRC': 0.0, '97LRD': 0.0, '97LRF': 0.0, '97LRG': 0.0, '97LRCRF01_AE': 0.0, '97LRCRF02_AG': 0.0, '97LROther': 0.0, '97LRAll': 0.0, '98NIA': 0.0, '98NIB': 0.4, '98NIC': 0.2, '98NID': 0.8, '98NIF': 0.1, '98NIG': 0.0, '98NICRF01_AE': 0.0, '98NICRF02_AG': 0.0, '98NIOther': 0.1, '98NIAll': 0.3, '98NKA': 0.0, '98NKB': 0.2, '98NKC': 0.1, '98NKD': 0.0, '98NKF': 0.1, '98NKG': 0.0, '98NKCRF01_AE': 0.0, '98NKCRF02_AG': 0.0, '98NKOther': 0.2, '98NKAll': 0.1, '98NDA': 0.3, '98NDB': 0.0, '98NDC': 0.0, '98NDD': 0.0, '98NDF': 0.0, '98NDG': 0.3, '98NDCRF01_AE': 0.0, '98NDCRF02_AG': 0.0, '98NDOther': 0.0, '98NDAll': 0.0, '98NSA': 0.0, '98NSB': 0.0, '98NSC': 0.0, '98NSD': 0.0, '98NSF': 0.0, '98NSG': 0.0, '98NSCRF01_AE': 0.2, '98NSCRF02_AG': 0.0, '98NSOther': 0.0, '98NSAll': 0.0, '98NHA': 0.0, '98NHB': 0.0, '98NHC': 0.0, '98NHD': 0.0, '98NHF': 0.0, '98NHG': 0.0, '98NHCRF01_AE': 0.0, '98NHCRF02_AG': 0.0, '98NHOther': 0.0, '98NHAll': 0.0, '98NYA': 0.0, '98NYB': 0.0, '98NYC': 0.0, '98NYD': 0.0, '98NYF': 0.0, '98NYG': 0.0, '98NYCRF01_AE': 0.0, '98NYCRF02_AG': 0.0, '98NYOther': 0.0, '98NYAll': 0.0, '98NFA': 0.0, '98NFB': 0.0, '98NFC': 0.0, '98NFD': 0.0, '98NFF': 0.0, '98NFG': 0.0, '98NFCRF01_AE': 0.0, '98NFCRF02_AG': 0.0, '98NFOther': 0.0, '98NFAll': 0.0, '98NTA': 0.0, '98NTB': 0.0, '98NTC': 0.0, '98NTD': 0.0, '98NTF': 0.0, '98NTG': 0.0, '98NTCRF01_AE': 0.0, '98NTCRF02_AG': 0.0, '98NTOther': 0.0, '98NTAll': 0.0, '98NQA': 0.0, '98NQB': 0.0, '98NQC': 0.0, '98NQD': 0.0, '98NQF': 0.0, '98NQG': 0.0, '98NQCRF01_AE': 0.0, '98NQCRF02_AG': 0.0, '98NQOther': 0.0, '98NQAll': 0.0, '98NCA': 0.0, '98NCB': 0.0, '98NCC': 0.0, '98NCD': 0.0, '98NCF': 0.0, '98NCG': 0.0, '98NCCRF01_AE': 0.0, '98NCCRF02_AG': 0.0, '98NCOther': 0.0, '98NCAll': 0.0, '98NRA': 0.0, '98NRB': 0.0, '98NRC': 0.0, '98NRD': 0.0, '98NRF': 0.0, '98NRG': 0.0, '98NRCRF01_AE': 0.0, '98NRCRF02_AG': 0.0, '98NROther': 0.0, '98NRAll': 0.0, '99FLA': 0.3, '99FLB': 0.1, '99FLC': 0.0, '99FLD': 0.0, '99FLF': 0.1, '99FLG': 0.8, '99FLCRF01_AE': 1.1, '99FLCRF02_AG': 0.4, '99FLOther': 0.2, '99FLAll': 0.1, '99FSA': 0.0, '99FSB': 0.0, '99FSC': 0.0, '99FSD': 0.0, '99FSF': 0.1, '99FSG': 0.0, '99FSCRF01_AE': 0.0, '99FSCRF02_AG': 0.0, '99FSOther': 0.0, '99FSAll': 0.0, '99FYA': 0.0, '99FYB': 0.1, '99FYC': 0.0, '99FYD': 0.0, '99FYF': 0.1, '99FYG': 0.0, '99FYCRF01_AE': 0.0, '99FYCRF02_AG': 0.0, '99FYOther': 0.0, '99FYAll': 0.0, '99FNA': 0.0, '99FNB': 0.0, '99FNC': 0.0, '99FND': 0.0, '99FNF': 0.1, '99FNG': 0.0, '99FNCRF01_AE': 0.0, '99FNCRF02_AG': 0.0, '99FNOther': 0.0, '99FNAll': 0.0, '99FVA': 0.0, '99FVB': 0.0, '99FVC': 0.0, '99FVD': 0.0, '99FVF': 0.0, '99FVG': 0.0, '99FVCRF01_AE': 0.0, '99FVCRF02_AG': 0.0, '99FVOther': 0.0, '99FVAll': 0.0, '99FCA': 0.0, '99FCB': 0.0, '99FCC': 0.0, '99FCD': 0.0, '99FCF': 0.0, '99FCG': 0.0, '99FCCRF01_AE': 0.0, '99FCCRF02_AG': 0.0, '99FCOther': 0.0, '99FCAll': 0.0, '99FWA': 0.0, '99FWB': 0.0, '99FWC': 0.0, '99FWD': 0.0, '99FWF': 0.0, '99FWG': 0.0, '99FWCRF01_AE': 0.0, '99FWCRF02_AG': 0.0, '99FWOther': 0.0, '99FWAll': 0.0, '99FIA': 0.0, '99FIB': 0.0, '99FIC': 0.0, '99FID': 0.0, '99FIF': 0.0, '99FIG': 0.0, '99FICRF01_AE': 0.0, '99FICRF02_AG': 0.0, '99FIOther': 0.0, '99FIAll': 0.0, '99FGA': 0.0, '99FGB': 0.0, '99FGC': 0.0, '99FGD': 0.0, '99FGF': 0.0, '99FGG': 0.0, '99FGCRF01_AE': 0.0, '99FGCRF02_AG': 0.0, '99FGOther': 0.0, '99FGAll': 0.0}
        result = self.aligner.prevalence_parser(filename)

        self.assertEqual(expected, result)

    def testGetAlignedSeq(self):
        # Setting params
        sites = [{'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
                 {'LengthNA': 3, 'PosAA': 59, 'PosNA': 7},
                 {'LengthNA': 3, 'PosAA': 60, 'PosNA': 10},
                 {'LengthNA': 3, 'PosAA': 61, 'PosNA': 13},
                 {'LengthNA': 3, 'PosAA': 62, 'PosNA': 16},
                 {'LengthNA': 3, 'PosAA': 63, 'PosNA': 19},
                 {'LengthNA': 3, 'PosAA': 64, 'PosNA': 22},
                 {'LengthNA': 3, 'PosAA': 65, 'PosNA': 25},
                 {'LengthNA': 3, 'PosAA': 66, 'PosNA': 28},
                 {'LengthNA': 3, 'PosAA': 67, 'PosNA': 31},
                 {'LengthNA': 3, 'PosAA': 68, 'PosNA': 34},
                 {'LengthNA': 3, 'PosAA': 69, 'PosNA': 37},
                 {'LengthNA': 3, 'PosAA': 70, 'PosNA': 40},
                 {'LengthNA': 3, 'PosAA': 71, 'PosNA': 43},
                 {'LengthNA': 3, 'PosAA': 72, 'PosNA': 46},
                 {'LengthNA': 3, 'PosAA': 73, 'PosNA': 49},
                 {'LengthNA': 3, 'PosAA': 74, 'PosNA': 52},
                 {'LengthNA': 3, 'PosAA': 75, 'PosNA': 55},
                 {'LengthNA': 3, 'PosAA': 76, 'PosNA': 58},
                 {'LengthNA': 3, 'PosAA': 77, 'PosNA': 61},
                 {'LengthNA': 3, 'PosAA': 78, 'PosNA': 64},
                 {'LengthNA': 3, 'PosAA': 79, 'PosNA': 67},
                 {'LengthNA': 3, 'PosAA': 80, 'PosNA': 70},
                 {'LengthNA': 3, 'PosAA': 81, 'PosNA': 73},
                 {'LengthNA': 3, 'PosAA': 82, 'PosNA': 76},
                 {'LengthNA': 3, 'PosAA': 83, 'PosNA': 79},
                 {'LengthNA': 3, 'PosAA': 84, 'PosNA': 82},
                 {'LengthNA': 3, 'PosAA': 85, 'PosNA': 85},
                 {'LengthNA': 3, 'PosAA': 86, 'PosNA': 88},
                 {'LengthNA': 3, 'PosAA': 87, 'PosNA': 91},
                 {'LengthNA': 3, 'PosAA': 88, 'PosNA': 94},
                 {'LengthNA': 3, 'PosAA': 89, 'PosNA': 97},
                 {'LengthNA': 3, 'PosAA': 90, 'PosNA': 100},
                 {'LengthNA': 3, 'PosAA': 91, 'PosNA': 103},
                 {'LengthNA': 3, 'PosAA': 92, 'PosNA': 106},
                 {'LengthNA': 3, 'PosAA': 93, 'PosNA': 109},
                 {'LengthNA': 3, 'PosAA': 94, 'PosNA': 112},
                 {'LengthNA': 3, 'PosAA': 95, 'PosNA': 115},
                 {'LengthNA': 3, 'PosAA': 96, 'PosNA': 118},
                 {'LengthNA': 3, 'PosAA': 97, 'PosNA': 121},
                 {'LengthNA': 3, 'PosAA': 98, 'PosNA': 124},
                 {'LengthNA': 3, 'PosAA': 99, 'PosNA': 127},
                 {'LengthNA': 3, 'PosAA': 100, 'PosNA': 130},
                 {'LengthNA': 3, 'PosAA': 101, 'PosNA': 133},
                 {'LengthNA': 3, 'PosAA': 102, 'PosNA': 136},
                 {'LengthNA': 3, 'PosAA': 103, 'PosNA': 139},
                 {'LengthNA': 3, 'PosAA': 104, 'PosNA': 142},
                 {'LengthNA': 3, 'PosAA': 105, 'PosNA': 145},
                 {'LengthNA': 3, 'PosAA': 106, 'PosNA': 148},
                 {'LengthNA': 3, 'PosAA': 107, 'PosNA': 151},
                 {'LengthNA': 3, 'PosAA': 108, 'PosNA': 154},
                 {'LengthNA': 3, 'PosAA': 109, 'PosNA': 157},
                 {'LengthNA': 3, 'PosAA': 110, 'PosNA': 160},
                 {'LengthNA': 3, 'PosAA': 111, 'PosNA': 163},
                 {'LengthNA': 3, 'PosAA': 112, 'PosNA': 166},
                 {'LengthNA': 3, 'PosAA': 113, 'PosNA': 169},
                 {'LengthNA': 3, 'PosAA': 114, 'PosNA': 172},
                 {'LengthNA': 3, 'PosAA': 115, 'PosNA': 175},
                 {'LengthNA': 3, 'PosAA': 116, 'PosNA': 178},
                 {'LengthNA': 3, 'PosAA': 117, 'PosNA': 181},
                 {'LengthNA': 3, 'PosAA': 118, 'PosNA': 184},
                 {'LengthNA': 3, 'PosAA': 119, 'PosNA': 187},
                 {'LengthNA': 3, 'PosAA': 120, 'PosNA': 190},
                 {'LengthNA': 3, 'PosAA': 121, 'PosNA': 193},
                 {'LengthNA': 3, 'PosAA': 122, 'PosNA': 196},
                 {'LengthNA': 3, 'PosAA': 123, 'PosNA': 199},
                 {'LengthNA': 3, 'PosAA': 124, 'PosNA': 202},
                 {'LengthNA': 3, 'PosAA': 125, 'PosNA': 205},
                 {'LengthNA': 3, 'PosAA': 126, 'PosNA': 208},
                 {'LengthNA': 3, 'PosAA': 127, 'PosNA': 211},
                 {'LengthNA': 3, 'PosAA': 128, 'PosNA': 214},
                 {'LengthNA': 3, 'PosAA': 129, 'PosNA': 217},
                 {'LengthNA': 3, 'PosAA': 130, 'PosNA': 220},
                 {'LengthNA': 3, 'PosAA': 131, 'PosNA': 223},
                 {'LengthNA': 3, 'PosAA': 132, 'PosNA': 226},
                 {'LengthNA': 3, 'PosAA': 133, 'PosNA': 229},
                 {'LengthNA': 3, 'PosAA': 134, 'PosNA': 232},
                 {'LengthNA': 3, 'PosAA': 135, 'PosNA': 235},
                 {'LengthNA': 3, 'PosAA': 136, 'PosNA': 238},
                 {'LengthNA': 3, 'PosAA': 137, 'PosNA': 241},
                 {'LengthNA': 3, 'PosAA': 138, 'PosNA': 244},
                 {'LengthNA': 3, 'PosAA': 139, 'PosNA': 247},
                 {'LengthNA': 3, 'PosAA': 140, 'PosNA': 250},
                 {'LengthNA': 3, 'PosAA': 141, 'PosNA': 253},
                 {'LengthNA': 3, 'PosAA': 142, 'PosNA': 256},
                 {'LengthNA': 3, 'PosAA': 143, 'PosNA': 259},
                 {'LengthNA': 3, 'PosAA': 144, 'PosNA': 262},
                 {'LengthNA': 3, 'PosAA': 145, 'PosNA': 265},
                 {'LengthNA': 3, 'PosAA': 146, 'PosNA': 268},
                 {'LengthNA': 3, 'PosAA': 147, 'PosNA': 271},
                 {'LengthNA': 3, 'PosAA': 148, 'PosNA': 274},
                 {'LengthNA': 3, 'PosAA': 149, 'PosNA': 277},
                 {'LengthNA': 3, 'PosAA': 150, 'PosNA': 280},
                 {'LengthNA': 3, 'PosAA': 151, 'PosNA': 283},
                 {'LengthNA': 3, 'PosAA': 152, 'PosNA': 286},
                 {'LengthNA': 3, 'PosAA': 153, 'PosNA': 289},
                 {'LengthNA': 3, 'PosAA': 154, 'PosNA': 292},
                 {'LengthNA': 3, 'PosAA': 155, 'PosNA': 295}]
        nuc = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'

        exp_aligned = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAAT'
        res_aligned = self.aligner.get_aligned_seq(nuc, sites)

        self.assertEqual(exp_aligned, res_aligned)

    def testGetAlignedSeqPost(self):
        filename = r'tests/hxb2-pr.fa'
        # Setting params
        sites = [{'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
                 {'LengthNA': 3, 'PosAA': 59, 'PosNA': 7},
                 {'LengthNA': 3, 'PosAA': 60, 'PosNA': 10},
                 {'LengthNA': 3, 'PosAA': 61, 'PosNA': 13},
                 {'LengthNA': 3, 'PosAA': 62, 'PosNA': 16},
                 {'LengthNA': 3, 'PosAA': 63, 'PosNA': 19},
                 {'LengthNA': 3, 'PosAA': 64, 'PosNA': 22},
                 {'LengthNA': 3, 'PosAA': 65, 'PosNA': 25},
                 {'LengthNA': 3, 'PosAA': 66, 'PosNA': 28},
                 {'LengthNA': 3, 'PosAA': 67, 'PosNA': 31},
                 {'LengthNA': 3, 'PosAA': 68, 'PosNA': 34},
                 {'LengthNA': 3, 'PosAA': 69, 'PosNA': 37},
                 {'LengthNA': 3, 'PosAA': 70, 'PosNA': 40},
                 {'LengthNA': 3, 'PosAA': 71, 'PosNA': 43},
                 {'LengthNA': 3, 'PosAA': 72, 'PosNA': 46},
                 {'LengthNA': 3, 'PosAA': 73, 'PosNA': 49},
                 {'LengthNA': 3, 'PosAA': 74, 'PosNA': 52},
                 {'LengthNA': 3, 'PosAA': 75, 'PosNA': 55},
                 {'LengthNA': 3, 'PosAA': 76, 'PosNA': 58},
                 {'LengthNA': 3, 'PosAA': 77, 'PosNA': 61},
                 {'LengthNA': 3, 'PosAA': 78, 'PosNA': 64},
                 {'LengthNA': 3, 'PosAA': 79, 'PosNA': 67},
                 {'LengthNA': 3, 'PosAA': 80, 'PosNA': 70},
                 {'LengthNA': 3, 'PosAA': 81, 'PosNA': 73},
                 {'LengthNA': 3, 'PosAA': 82, 'PosNA': 76},
                 {'LengthNA': 3, 'PosAA': 83, 'PosNA': 79},
                 {'LengthNA': 3, 'PosAA': 84, 'PosNA': 82},
                 {'LengthNA': 3, 'PosAA': 85, 'PosNA': 85},
                 {'LengthNA': 3, 'PosAA': 86, 'PosNA': 88},
                 {'LengthNA': 3, 'PosAA': 87, 'PosNA': 91},
                 {'LengthNA': 3, 'PosAA': 88, 'PosNA': 94},
                 {'LengthNA': 3, 'PosAA': 89, 'PosNA': 97},
                 {'LengthNA': 3, 'PosAA': 90, 'PosNA': 100},
                 {'LengthNA': 3, 'PosAA': 91, 'PosNA': 103},
                 {'LengthNA': 3, 'PosAA': 92, 'PosNA': 106},
                 {'LengthNA': 3, 'PosAA': 93, 'PosNA': 109},
                 {'LengthNA': 3, 'PosAA': 94, 'PosNA': 112},
                 {'LengthNA': 3, 'PosAA': 95, 'PosNA': 115},
                 {'LengthNA': 3, 'PosAA': 96, 'PosNA': 118},
                 {'LengthNA': 3, 'PosAA': 97, 'PosNA': 121},
                 {'LengthNA': 3, 'PosAA': 98, 'PosNA': 124},
                 {'LengthNA': 3, 'PosAA': 99, 'PosNA': 127},
                 {'LengthNA': 3, 'PosAA': 100, 'PosNA': 130},
                 {'LengthNA': 3, 'PosAA': 101, 'PosNA': 133},
                 {'LengthNA': 3, 'PosAA': 102, 'PosNA': 136},
                 {'LengthNA': 3, 'PosAA': 103, 'PosNA': 139},
                 {'LengthNA': 3, 'PosAA': 104, 'PosNA': 142},
                 {'LengthNA': 3, 'PosAA': 105, 'PosNA': 145},
                 {'LengthNA': 3, 'PosAA': 106, 'PosNA': 148},
                 {'LengthNA': 3, 'PosAA': 107, 'PosNA': 151},
                 {'LengthNA': 3, 'PosAA': 108, 'PosNA': 154},
                 {'LengthNA': 3, 'PosAA': 109, 'PosNA': 157},
                 {'LengthNA': 3, 'PosAA': 110, 'PosNA': 160},
                 {'LengthNA': 3, 'PosAA': 111, 'PosNA': 163},
                 {'LengthNA': 3, 'PosAA': 112, 'PosNA': 166},
                 {'LengthNA': 3, 'PosAA': 113, 'PosNA': 169},
                 {'LengthNA': 3, 'PosAA': 114, 'PosNA': 172},
                 {'LengthNA': 3, 'PosAA': 115, 'PosNA': 175},
                 {'LengthNA': 3, 'PosAA': 116, 'PosNA': 178},
                 {'LengthNA': 3, 'PosAA': 117, 'PosNA': 181},
                 {'LengthNA': 3, 'PosAA': 118, 'PosNA': 184},
                 {'LengthNA': 3, 'PosAA': 119, 'PosNA': 187},
                 {'LengthNA': 3, 'PosAA': 120, 'PosNA': 190},
                 {'LengthNA': 3, 'PosAA': 121, 'PosNA': 193},
                 {'LengthNA': 3, 'PosAA': 122, 'PosNA': 196},
                 {'LengthNA': 3, 'PosAA': 123, 'PosNA': 199},
                 {'LengthNA': 3, 'PosAA': 124, 'PosNA': 202},
                 {'LengthNA': 3, 'PosAA': 125, 'PosNA': 205},
                 {'LengthNA': 3, 'PosAA': 126, 'PosNA': 208},
                 {'LengthNA': 3, 'PosAA': 127, 'PosNA': 211},
                 {'LengthNA': 3, 'PosAA': 128, 'PosNA': 214},
                 {'LengthNA': 3, 'PosAA': 129, 'PosNA': 217},
                 {'LengthNA': 3, 'PosAA': 130, 'PosNA': 220},
                 {'LengthNA': 3, 'PosAA': 131, 'PosNA': 223},
                 {'LengthNA': 3, 'PosAA': 132, 'PosNA': 226},
                 {'LengthNA': 3, 'PosAA': 133, 'PosNA': 229},
                 {'LengthNA': 3, 'PosAA': 134, 'PosNA': 232},
                 {'LengthNA': 3, 'PosAA': 135, 'PosNA': 235},
                 {'LengthNA': 3, 'PosAA': 136, 'PosNA': 238},
                 {'LengthNA': 3, 'PosAA': 137, 'PosNA': 241},
                 {'LengthNA': 3, 'PosAA': 138, 'PosNA': 244},
                 {'LengthNA': 3, 'PosAA': 139, 'PosNA': 247},
                 {'LengthNA': 3, 'PosAA': 140, 'PosNA': 250},
                 {'LengthNA': 3, 'PosAA': 141, 'PosNA': 253},
                 {'LengthNA': 3, 'PosAA': 142, 'PosNA': 256},
                 {'LengthNA': 3, 'PosAA': 143, 'PosNA': 259},
                 {'LengthNA': 3, 'PosAA': 144, 'PosNA': 262},
                 {'LengthNA': 3, 'PosAA': 145, 'PosNA': 265},
                 {'LengthNA': 3, 'PosAA': 146, 'PosNA': 268},
                 {'LengthNA': 3, 'PosAA': 147, 'PosNA': 271},
                 {'LengthNA': 3, 'PosAA': 148, 'PosNA': 274},
                 {'LengthNA': 3, 'PosAA': 149, 'PosNA': 277},
                 {'LengthNA': 3, 'PosAA': 150, 'PosNA': 280},
                 {'LengthNA': 3, 'PosAA': 151, 'PosNA': 283},
                 {'LengthNA': 3, 'PosAA': 152, 'PosNA': 286},
                 {'LengthNA': 3, 'PosAA': 153, 'PosNA': 289},
                 {'LengthNA': 3, 'PosAA': 154, 'PosNA': 292},
                 {'LengthNA': 3, 'PosAA': 155, 'PosNA': 295}]
        nuc = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'

        exp_aligned = 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAAT'
        res_aligned = self.postAligner.get_aligned_seq(nuc, sites)

        self.assertEqual(exp_aligned, res_aligned)

    def testAlignFile(self):
        # Setting params
        filename = r'tests/hxb2-pr.fa'

        exp_records = \
            [{'AlignedSites': [{'LengthNA': 3, 'PosAA': 57, 'PosNA': 1},
                               {'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
                               {'LengthNA': 3, 'PosAA': 59, 'PosNA': 7},
                               {'LengthNA': 3, 'PosAA': 60, 'PosNA': 10},
                               {'LengthNA': 3, 'PosAA': 61, 'PosNA': 13},
                               {'LengthNA': 3, 'PosAA': 62, 'PosNA': 16},
                               {'LengthNA': 3, 'PosAA': 63, 'PosNA': 19},
                               {'LengthNA': 3, 'PosAA': 64, 'PosNA': 22},
                               {'LengthNA': 3, 'PosAA': 65, 'PosNA': 25},
                               {'LengthNA': 3, 'PosAA': 66, 'PosNA': 28},
                               {'LengthNA': 3, 'PosAA': 67, 'PosNA': 31},
                               {'LengthNA': 3, 'PosAA': 68, 'PosNA': 34},
                               {'LengthNA': 3, 'PosAA': 69, 'PosNA': 37},
                               {'LengthNA': 3, 'PosAA': 70, 'PosNA': 40},
                               {'LengthNA': 3, 'PosAA': 71, 'PosNA': 43},
                               {'LengthNA': 3, 'PosAA': 72, 'PosNA': 46},
                               {'LengthNA': 3, 'PosAA': 73, 'PosNA': 49},
                               {'LengthNA': 3, 'PosAA': 74, 'PosNA': 52},
                               {'LengthNA': 3, 'PosAA': 75, 'PosNA': 55},
                               {'LengthNA': 3, 'PosAA': 76, 'PosNA': 58},
                               {'LengthNA': 3, 'PosAA': 77, 'PosNA': 61},
                               {'LengthNA': 3, 'PosAA': 78, 'PosNA': 64},
                               {'LengthNA': 3, 'PosAA': 79, 'PosNA': 67},
                               {'LengthNA': 3, 'PosAA': 80, 'PosNA': 70},
                               {'LengthNA': 3, 'PosAA': 81, 'PosNA': 73},
                               {'LengthNA': 3, 'PosAA': 82, 'PosNA': 76},
                               {'LengthNA': 3, 'PosAA': 83, 'PosNA': 79},
                               {'LengthNA': 3, 'PosAA': 84, 'PosNA': 82},
                               {'LengthNA': 3, 'PosAA': 85, 'PosNA': 85},
                               {'LengthNA': 3, 'PosAA': 86, 'PosNA': 88},
                               {'LengthNA': 3, 'PosAA': 87, 'PosNA': 91},
                               {'LengthNA': 3, 'PosAA': 88, 'PosNA': 94},
                               {'LengthNA': 3, 'PosAA': 89, 'PosNA': 97},
                               {'LengthNA': 3, 'PosAA': 90, 'PosNA': 100},
                               {'LengthNA': 3, 'PosAA': 91, 'PosNA': 103},
                               {'LengthNA': 3, 'PosAA': 92, 'PosNA': 106},
                               {'LengthNA': 3, 'PosAA': 93, 'PosNA': 109},
                               {'LengthNA': 3, 'PosAA': 94, 'PosNA': 112},
                               {'LengthNA': 3, 'PosAA': 95, 'PosNA': 115},
                               {'LengthNA': 3, 'PosAA': 96, 'PosNA': 118},
                               {'LengthNA': 3, 'PosAA': 97, 'PosNA': 121},
                               {'LengthNA': 3, 'PosAA': 98, 'PosNA': 124},
                               {'LengthNA': 3, 'PosAA': 99, 'PosNA': 127},
                               {'LengthNA': 3, 'PosAA': 100, 'PosNA': 130},
                               {'LengthNA': 3, 'PosAA': 101, 'PosNA': 133},
                               {'LengthNA': 3, 'PosAA': 102, 'PosNA': 136},
                               {'LengthNA': 3, 'PosAA': 103, 'PosNA': 139},
                               {'LengthNA': 3, 'PosAA': 104, 'PosNA': 142},
                               {'LengthNA': 3, 'PosAA': 105, 'PosNA': 145},
                               {'LengthNA': 3, 'PosAA': 106, 'PosNA': 148},
                               {'LengthNA': 3, 'PosAA': 107, 'PosNA': 151},
                               {'LengthNA': 3, 'PosAA': 108, 'PosNA': 154},
                               {'LengthNA': 3, 'PosAA': 109, 'PosNA': 157},
                               {'LengthNA': 3, 'PosAA': 110, 'PosNA': 160},
                               {'LengthNA': 3, 'PosAA': 111, 'PosNA': 163},
                               {'LengthNA': 3, 'PosAA': 112, 'PosNA': 166},
                               {'LengthNA': 3, 'PosAA': 113, 'PosNA': 169},
                               {'LengthNA': 3, 'PosAA': 114, 'PosNA': 172},
                               {'LengthNA': 3, 'PosAA': 115, 'PosNA': 175},
                               {'LengthNA': 3, 'PosAA': 116, 'PosNA': 178},
                               {'LengthNA': 3, 'PosAA': 117, 'PosNA': 181},
                               {'LengthNA': 3, 'PosAA': 118, 'PosNA': 184},
                               {'LengthNA': 3, 'PosAA': 119, 'PosNA': 187},
                               {'LengthNA': 3, 'PosAA': 120, 'PosNA': 190},
                               {'LengthNA': 3, 'PosAA': 121, 'PosNA': 193},
                               {'LengthNA': 3, 'PosAA': 122, 'PosNA': 196},
                               {'LengthNA': 3, 'PosAA': 123, 'PosNA': 199},
                               {'LengthNA': 3, 'PosAA': 124, 'PosNA': 202},
                               {'LengthNA': 3, 'PosAA': 125, 'PosNA': 205},
                               {'LengthNA': 3, 'PosAA': 126, 'PosNA': 208},
                               {'LengthNA': 3, 'PosAA': 127, 'PosNA': 211},
                               {'LengthNA': 3, 'PosAA': 128, 'PosNA': 214},
                               {'LengthNA': 3, 'PosAA': 129, 'PosNA': 217},
                               {'LengthNA': 3, 'PosAA': 130, 'PosNA': 220},
                               {'LengthNA': 3, 'PosAA': 131, 'PosNA': 223},
                               {'LengthNA': 3, 'PosAA': 132, 'PosNA': 226},
                               {'LengthNA': 3, 'PosAA': 133, 'PosNA': 229},
                               {'LengthNA': 3, 'PosAA': 134, 'PosNA': 232},
                               {'LengthNA': 3, 'PosAA': 135, 'PosNA': 235},
                               {'LengthNA': 3, 'PosAA': 136, 'PosNA': 238},
                               {'LengthNA': 3, 'PosAA': 137, 'PosNA': 241},
                               {'LengthNA': 3, 'PosAA': 138, 'PosNA': 244},
                               {'LengthNA': 3, 'PosAA': 139, 'PosNA': 247},
                               {'LengthNA': 3, 'PosAA': 140, 'PosNA': 250},
                               {'LengthNA': 3, 'PosAA': 141, 'PosNA': 253},
                               {'LengthNA': 3, 'PosAA': 142, 'PosNA': 256},
                               {'LengthNA': 3, 'PosAA': 143, 'PosNA': 259},
                               {'LengthNA': 3, 'PosAA': 144, 'PosNA': 262},
                               {'LengthNA': 3, 'PosAA': 145, 'PosNA': 265},
                               {'LengthNA': 3, 'PosAA': 146, 'PosNA': 268},
                               {'LengthNA': 3, 'PosAA': 147, 'PosNA': 271},
                               {'LengthNA': 3, 'PosAA': 148, 'PosNA': 274},
                               {'LengthNA': 3, 'PosAA': 149, 'PosNA': 277},
                               {'LengthNA': 3, 'PosAA': 150, 'PosNA': 280},
                               {'LengthNA': 3, 'PosAA': 151, 'PosNA': 283},
                               {'LengthNA': 3, 'PosAA': 152, 'PosNA': 286},
                               {'LengthNA': 3, 'PosAA': 153, 'PosNA': 289},
                               {'LengthNA': 3, 'PosAA': 154, 'PosNA': 292},
                               {'LengthNA': 3, 'PosAA': 155, 'PosNA': 295}],
              'FirstAA': 57,
              'FirstNA': 1,
              'Frameshifts': [],
              'LastAA': 155,
              'LastNA': 297,
              'Mutations': [{'AminoAcidText': 'V',
                             'CodonText': 'GTC',
                             'Control': '...',
                             'InsertedAminoAcidsText': '',
                             'InsertedCodonsText': '',
                             'IsDeletion': False,
                             'IsInsertion': False,
                             'IsPartial': False,
                             'NAPosition': 7,
                             'Position': 59,
                             'ReferenceText': 'I'},
                            {'AminoAcidText': 'S',
                             'CodonText': 'AGT',
                             'Control': '...',
                             'InsertedAminoAcidsText': '',
                             'InsertedCodonsText': '',
                             'IsDeletion': False,
                             'IsInsertion': False,
                             'IsPartial': False,
                             'NAPosition': 109,
                             'Position': 93,
                             'ReferenceText': 'N'}],
             'Name': 'HXB2-PR',
             'Sequence': 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
            {'AlignedSites': [{'LengthNA': 3, 'PosAA': 58, 'PosNA': 3},
                              {'LengthNA': 3, 'PosAA': 59, 'PosNA': 6},
                              {'LengthNA': 3, 'PosAA': 60, 'PosNA': 9},
                              {'LengthNA': 3, 'PosAA': 61, 'PosNA': 12},
                              {'LengthNA': 3, 'PosAA': 62, 'PosNA': 15},
                              {'LengthNA': 3, 'PosAA': 63, 'PosNA': 18},
                              {'LengthNA': 3, 'PosAA': 64, 'PosNA': 21},
                              {'LengthNA': 3, 'PosAA': 65, 'PosNA': 24},
                              {'LengthNA': 3, 'PosAA': 66, 'PosNA': 27},
                              {'LengthNA': 3, 'PosAA': 67, 'PosNA': 30},
                              {'LengthNA': 3, 'PosAA': 68, 'PosNA': 33},
                              {'LengthNA': 3, 'PosAA': 69, 'PosNA': 36},
                              {'LengthNA': 3, 'PosAA': 70, 'PosNA': 39},
                              {'LengthNA': 3, 'PosAA': 71, 'PosNA': 42},
                              {'LengthNA': 3, 'PosAA': 72, 'PosNA': 45},
                              {'LengthNA': 3, 'PosAA': 73, 'PosNA': 48},
                              {'LengthNA': 3, 'PosAA': 74, 'PosNA': 51},
                              {'LengthNA': 3, 'PosAA': 75, 'PosNA': 54},
                              {'LengthNA': 3, 'PosAA': 76, 'PosNA': 57},
                              {'LengthNA': 3, 'PosAA': 77, 'PosNA': 60},
                              {'LengthNA': 3, 'PosAA': 78, 'PosNA': 63},
                              {'LengthNA': 3, 'PosAA': 79, 'PosNA': 66},
                              {'LengthNA': 3, 'PosAA': 80, 'PosNA': 69},
                              {'LengthNA': 3, 'PosAA': 81, 'PosNA': 72},
                              {'LengthNA': 3, 'PosAA': 82, 'PosNA': 75},
                              {'LengthNA': 3, 'PosAA': 83, 'PosNA': 78},
                              {'LengthNA': 3, 'PosAA': 84, 'PosNA': 81},
                              {'LengthNA': 3, 'PosAA': 85, 'PosNA': 84},
                              {'LengthNA': 3, 'PosAA': 86, 'PosNA': 87},
                              {'LengthNA': 3, 'PosAA': 87, 'PosNA': 90},
                              {'LengthNA': 3, 'PosAA': 88, 'PosNA': 93},
                              {'LengthNA': 3, 'PosAA': 89, 'PosNA': 96},
                              {'LengthNA': 3, 'PosAA': 90, 'PosNA': 99},
                              {'LengthNA': 3, 'PosAA': 91, 'PosNA': 102},
                              {'LengthNA': 3, 'PosAA': 92, 'PosNA': 105},
                              {'LengthNA': 3, 'PosAA': 93, 'PosNA': 108},
                              {'LengthNA': 3, 'PosAA': 94, 'PosNA': 111},
                              {'LengthNA': 3, 'PosAA': 95, 'PosNA': 114},
                              {'LengthNA': 3, 'PosAA': 96, 'PosNA': 117},
                              {'LengthNA': 3, 'PosAA': 97, 'PosNA': 120},
                              {'LengthNA': 3, 'PosAA': 98, 'PosNA': 123},
                              {'LengthNA': 3, 'PosAA': 99, 'PosNA': 126},
                              {'LengthNA': 3, 'PosAA': 100, 'PosNA': 129},
                              {'LengthNA': 3, 'PosAA': 101, 'PosNA': 132},
                              {'LengthNA': 3, 'PosAA': 102, 'PosNA': 135},
                              {'LengthNA': 3, 'PosAA': 103, 'PosNA': 138},
                              {'LengthNA': 3, 'PosAA': 104, 'PosNA': 141},
                              {'LengthNA': 3, 'PosAA': 105, 'PosNA': 144},
                              {'LengthNA': 3, 'PosAA': 106, 'PosNA': 147},
                              {'LengthNA': 3, 'PosAA': 107, 'PosNA': 150},
                              {'LengthNA': 3, 'PosAA': 108, 'PosNA': 153},
                              {'LengthNA': 3, 'PosAA': 109, 'PosNA': 156},
                              {'LengthNA': 3, 'PosAA': 110, 'PosNA': 159},
                              {'LengthNA': 3, 'PosAA': 111, 'PosNA': 162},
                              {'LengthNA': 3, 'PosAA': 112, 'PosNA': 165},
                              {'LengthNA': 3, 'PosAA': 113, 'PosNA': 168},
                              {'LengthNA': 3, 'PosAA': 114, 'PosNA': 171},
                              {'LengthNA': 3, 'PosAA': 115, 'PosNA': 174},
                              {'LengthNA': 3, 'PosAA': 116, 'PosNA': 177},
                              {'LengthNA': 3, 'PosAA': 117, 'PosNA': 180},
                              {'LengthNA': 3, 'PosAA': 118, 'PosNA': 183},
                              {'LengthNA': 3, 'PosAA': 119, 'PosNA': 186},
                              {'LengthNA': 3, 'PosAA': 120, 'PosNA': 189},
                              {'LengthNA': 3, 'PosAA': 121, 'PosNA': 192},
                              {'LengthNA': 3, 'PosAA': 122, 'PosNA': 195},
                              {'LengthNA': 3, 'PosAA': 123, 'PosNA': 198},
                              {'LengthNA': 3, 'PosAA': 124, 'PosNA': 201},
                              {'LengthNA': 3, 'PosAA': 125, 'PosNA': 204},
                              {'LengthNA': 3, 'PosAA': 126, 'PosNA': 207},
                              {'LengthNA': 3, 'PosAA': 127, 'PosNA': 210},
                              {'LengthNA': 3, 'PosAA': 128, 'PosNA': 213},
                              {'LengthNA': 3, 'PosAA': 129, 'PosNA': 216},
                              {'LengthNA': 3, 'PosAA': 130, 'PosNA': 219},
                              {'LengthNA': 3, 'PosAA': 131, 'PosNA': 222},
                              {'LengthNA': 3, 'PosAA': 132, 'PosNA': 225},
                              {'LengthNA': 3, 'PosAA': 133, 'PosNA': 228},
                              {'LengthNA': 3, 'PosAA': 134, 'PosNA': 231},
                              {'LengthNA': 3, 'PosAA': 135, 'PosNA': 234},
                              {'LengthNA': 3, 'PosAA': 136, 'PosNA': 237},
                              {'LengthNA': 3, 'PosAA': 137, 'PosNA': 240},
                              {'LengthNA': 3, 'PosAA': 138, 'PosNA': 243},
                              {'LengthNA': 3, 'PosAA': 139, 'PosNA': 246},
                              {'LengthNA': 3, 'PosAA': 140, 'PosNA': 249},
                              {'LengthNA': 3, 'PosAA': 141, 'PosNA': 252},
                              {'LengthNA': 3, 'PosAA': 142, 'PosNA': 255},
                              {'LengthNA': 3, 'PosAA': 143, 'PosNA': 258},
                              {'LengthNA': 3, 'PosAA': 144, 'PosNA': 261},
                              {'LengthNA': 3, 'PosAA': 145, 'PosNA': 264},
                              {'LengthNA': 3, 'PosAA': 146, 'PosNA': 267},
                              {'LengthNA': 3, 'PosAA': 147, 'PosNA': 270},
                              {'LengthNA': 3, 'PosAA': 148, 'PosNA': 273},
                              {'LengthNA': 3, 'PosAA': 149, 'PosNA': 276},
                              {'LengthNA': 3, 'PosAA': 150, 'PosNA': 279},
                              {'LengthNA': 3, 'PosAA': 151, 'PosNA': 282},
                              {'LengthNA': 3, 'PosAA': 152, 'PosNA': 285},
                              {'LengthNA': 3, 'PosAA': 153, 'PosNA': 288},
                              {'LengthNA': 3, 'PosAA': 154, 'PosNA': 291},
                              {'LengthNA': 3, 'PosAA': 155, 'PosNA': 294}],
            'FirstAA': 58,
            'FirstNA': 3,
            'Frameshifts': [],
            'LastAA': 155,
            'LastNA': 296,
            'Mutations': [{'AminoAcidText': 'V',
                            'CodonText': 'GTC',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 6,
                            'Position': 59,
                            'ReferenceText': 'I'},
                            {'AminoAcidText': 'S',
                            'CodonText': 'AGT',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 108,
                            'Position': 93,
                            'ReferenceText': 'N'}],
            'Name': 'shift1',
            'Sequence': 'CAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
            {'AlignedSites': [{'LengthNA': 3, 'PosAA': 58, 'PosNA': 2},
                              {'LengthNA': 3, 'PosAA': 59, 'PosNA': 5},
                              {'LengthNA': 3, 'PosAA': 60, 'PosNA': 8},
                              {'LengthNA': 3, 'PosAA': 61, 'PosNA': 11},
                              {'LengthNA': 3, 'PosAA': 62, 'PosNA': 14},
                              {'LengthNA': 3, 'PosAA': 63, 'PosNA': 17},
                              {'LengthNA': 3, 'PosAA': 64, 'PosNA': 20},
                              {'LengthNA': 3, 'PosAA': 65, 'PosNA': 23},
                              {'LengthNA': 3, 'PosAA': 66, 'PosNA': 26},
                              {'LengthNA': 3, 'PosAA': 67, 'PosNA': 29},
                              {'LengthNA': 3, 'PosAA': 68, 'PosNA': 32},
                              {'LengthNA': 3, 'PosAA': 69, 'PosNA': 35},
                              {'LengthNA': 3, 'PosAA': 70, 'PosNA': 38},
                              {'LengthNA': 3, 'PosAA': 71, 'PosNA': 41},
                              {'LengthNA': 3, 'PosAA': 72, 'PosNA': 44},
                              {'LengthNA': 3, 'PosAA': 73, 'PosNA': 47},
                              {'LengthNA': 3, 'PosAA': 74, 'PosNA': 50},
                              {'LengthNA': 3, 'PosAA': 75, 'PosNA': 53},
                              {'LengthNA': 3, 'PosAA': 76, 'PosNA': 56},
                              {'LengthNA': 3, 'PosAA': 77, 'PosNA': 59},
                              {'LengthNA': 3, 'PosAA': 78, 'PosNA': 62},
                              {'LengthNA': 3, 'PosAA': 79, 'PosNA': 65},
                              {'LengthNA': 3, 'PosAA': 80, 'PosNA': 68},
                              {'LengthNA': 3, 'PosAA': 81, 'PosNA': 71},
                              {'LengthNA': 3, 'PosAA': 82, 'PosNA': 74},
                              {'LengthNA': 3, 'PosAA': 83, 'PosNA': 77},
                              {'LengthNA': 3, 'PosAA': 84, 'PosNA': 80},
                              {'LengthNA': 3, 'PosAA': 85, 'PosNA': 83},
                              {'LengthNA': 3, 'PosAA': 86, 'PosNA': 86},
                              {'LengthNA': 3, 'PosAA': 87, 'PosNA': 89},
                              {'LengthNA': 3, 'PosAA': 88, 'PosNA': 92},
                              {'LengthNA': 3, 'PosAA': 89, 'PosNA': 95},
                              {'LengthNA': 3, 'PosAA': 90, 'PosNA': 98},
                              {'LengthNA': 3, 'PosAA': 91, 'PosNA': 101},
                              {'LengthNA': 3, 'PosAA': 92, 'PosNA': 104},
                              {'LengthNA': 3, 'PosAA': 93, 'PosNA': 107},
                              {'LengthNA': 3, 'PosAA': 94, 'PosNA': 110},
                              {'LengthNA': 3, 'PosAA': 95, 'PosNA': 113},
                              {'LengthNA': 3, 'PosAA': 96, 'PosNA': 116},
                              {'LengthNA': 3, 'PosAA': 97, 'PosNA': 119},
                              {'LengthNA': 3, 'PosAA': 98, 'PosNA': 122},
                              {'LengthNA': 3, 'PosAA': 99, 'PosNA': 125},
                              {'LengthNA': 3, 'PosAA': 100, 'PosNA': 128},
                              {'LengthNA': 3, 'PosAA': 101, 'PosNA': 131},
                              {'LengthNA': 3, 'PosAA': 102, 'PosNA': 134},
                              {'LengthNA': 3, 'PosAA': 103, 'PosNA': 137},
                              {'LengthNA': 3, 'PosAA': 104, 'PosNA': 140},
                              {'LengthNA': 3, 'PosAA': 105, 'PosNA': 143},
                              {'LengthNA': 3, 'PosAA': 106, 'PosNA': 146},
                              {'LengthNA': 3, 'PosAA': 107, 'PosNA': 149},
                              {'LengthNA': 3, 'PosAA': 108, 'PosNA': 152},
                              {'LengthNA': 3, 'PosAA': 109, 'PosNA': 155},
                              {'LengthNA': 3, 'PosAA': 110, 'PosNA': 158},
                              {'LengthNA': 3, 'PosAA': 111, 'PosNA': 161},
                              {'LengthNA': 3, 'PosAA': 112, 'PosNA': 164},
                              {'LengthNA': 3, 'PosAA': 113, 'PosNA': 167},
                              {'LengthNA': 3, 'PosAA': 114, 'PosNA': 170},
                              {'LengthNA': 3, 'PosAA': 115, 'PosNA': 173},
                              {'LengthNA': 3, 'PosAA': 116, 'PosNA': 176},
                              {'LengthNA': 3, 'PosAA': 117, 'PosNA': 179},
                              {'LengthNA': 3, 'PosAA': 118, 'PosNA': 182},
                              {'LengthNA': 3, 'PosAA': 119, 'PosNA': 185},
                              {'LengthNA': 3, 'PosAA': 120, 'PosNA': 188},
                              {'LengthNA': 3, 'PosAA': 121, 'PosNA': 191},
                              {'LengthNA': 3, 'PosAA': 122, 'PosNA': 194},
                              {'LengthNA': 3, 'PosAA': 123, 'PosNA': 197},
                              {'LengthNA': 3, 'PosAA': 124, 'PosNA': 200},
                              {'LengthNA': 3, 'PosAA': 125, 'PosNA': 203},
                              {'LengthNA': 3, 'PosAA': 126, 'PosNA': 206},
                              {'LengthNA': 3, 'PosAA': 127, 'PosNA': 209},
                              {'LengthNA': 3, 'PosAA': 128, 'PosNA': 212},
                              {'LengthNA': 3, 'PosAA': 129, 'PosNA': 215},
                              {'LengthNA': 3, 'PosAA': 130, 'PosNA': 218},
                              {'LengthNA': 3, 'PosAA': 131, 'PosNA': 221},
                              {'LengthNA': 3, 'PosAA': 132, 'PosNA': 224},
                              {'LengthNA': 3, 'PosAA': 133, 'PosNA': 227},
                              {'LengthNA': 3, 'PosAA': 134, 'PosNA': 230},
                              {'LengthNA': 3, 'PosAA': 135, 'PosNA': 233},
                              {'LengthNA': 3, 'PosAA': 136, 'PosNA': 236},
                              {'LengthNA': 3, 'PosAA': 137, 'PosNA': 239},
                              {'LengthNA': 3, 'PosAA': 138, 'PosNA': 242},
                              {'LengthNA': 3, 'PosAA': 139, 'PosNA': 245},
                              {'LengthNA': 3, 'PosAA': 140, 'PosNA': 248},
                              {'LengthNA': 3, 'PosAA': 141, 'PosNA': 251},
                              {'LengthNA': 3, 'PosAA': 142, 'PosNA': 254},
                              {'LengthNA': 3, 'PosAA': 143, 'PosNA': 257},
                              {'LengthNA': 3, 'PosAA': 144, 'PosNA': 260},
                              {'LengthNA': 3, 'PosAA': 145, 'PosNA': 263},
                              {'LengthNA': 3, 'PosAA': 146, 'PosNA': 266},
                              {'LengthNA': 3, 'PosAA': 147, 'PosNA': 269},
                              {'LengthNA': 3, 'PosAA': 148, 'PosNA': 272},
                              {'LengthNA': 3, 'PosAA': 149, 'PosNA': 275},
                              {'LengthNA': 3, 'PosAA': 150, 'PosNA': 278},
                              {'LengthNA': 3, 'PosAA': 151, 'PosNA': 281},
                              {'LengthNA': 3, 'PosAA': 152, 'PosNA': 284},
                              {'LengthNA': 3, 'PosAA': 153, 'PosNA': 287},
                              {'LengthNA': 3, 'PosAA': 154, 'PosNA': 290},
                              {'LengthNA': 3, 'PosAA': 155, 'PosNA': 293}],
              'FirstAA': 58,
              'FirstNA': 2,
              'Frameshifts': [],
              'LastAA': 155,
              'LastNA': 295,
              'Mutations': [{'AminoAcidText': 'V',
                             'CodonText': 'GTC',
                             'Control': '...',
                             'InsertedAminoAcidsText': '',
                             'InsertedCodonsText': '',
                             'IsDeletion': False,
                             'IsInsertion': False,
                             'IsPartial': False,
                             'NAPosition': 5,
                             'Position': 59,
                             'ReferenceText': 'I'},
                            {'AminoAcidText': 'S',
                             'CodonText': 'AGT',
                             'Control': '...',
                             'InsertedAminoAcidsText': '',
                             'InsertedCodonsText': '',
                             'IsDeletion': False,
                             'IsInsertion': False,
                             'IsPartial': False,
                             'NAPosition': 107,
                             'Position': 93,
                             'ReferenceText': 'N'}],
              'Name': 'shift2',
              'Sequence': 'CAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
             {'AlignedSites': [{'LengthNA': 3, 'PosAA': 57, 'PosNA': 2},
                               {'LengthNA': 3, 'PosAA': 58, 'PosNA': 5},
                               {'LengthNA': 3, 'PosAA': 59, 'PosNA': 8},
                               {'LengthNA': 3, 'PosAA': 60, 'PosNA': 11},
                               {'LengthNA': 3, 'PosAA': 61, 'PosNA': 14},
                               {'LengthNA': 3, 'PosAA': 62, 'PosNA': 17},
                               {'LengthNA': 3, 'PosAA': 63, 'PosNA': 20},
                               {'LengthNA': 3, 'PosAA': 64, 'PosNA': 23},
                               {'LengthNA': 3, 'PosAA': 65, 'PosNA': 26},
                               {'LengthNA': 3, 'PosAA': 66, 'PosNA': 29},
                               {'LengthNA': 3, 'PosAA': 67, 'PosNA': 32},
                               {'LengthNA': 3, 'PosAA': 68, 'PosNA': 35},
                               {'LengthNA': 3, 'PosAA': 69, 'PosNA': 38},
                               {'LengthNA': 3, 'PosAA': 70, 'PosNA': 41},
                               {'LengthNA': 3, 'PosAA': 71, 'PosNA': 44},
                               {'LengthNA': 3, 'PosAA': 72, 'PosNA': 47},
                               {'LengthNA': 3, 'PosAA': 73, 'PosNA': 50},
                               {'LengthNA': 3, 'PosAA': 74, 'PosNA': 53},
                               {'LengthNA': 3, 'PosAA': 75, 'PosNA': 56},
                               {'LengthNA': 3, 'PosAA': 76, 'PosNA': 59},
                               {'LengthNA': 3, 'PosAA': 77, 'PosNA': 62},
                               {'LengthNA': 3, 'PosAA': 78, 'PosNA': 65},
                               {'LengthNA': 3, 'PosAA': 79, 'PosNA': 68},
                               {'LengthNA': 3, 'PosAA': 80, 'PosNA': 71},
                               {'LengthNA': 3, 'PosAA': 81, 'PosNA': 74},
                               {'LengthNA': 3, 'PosAA': 82, 'PosNA': 77},
                               {'LengthNA': 3, 'PosAA': 83, 'PosNA': 80},
                               {'LengthNA': 3, 'PosAA': 84, 'PosNA': 83},
                               {'LengthNA': 3, 'PosAA': 85, 'PosNA': 86},
                               {'LengthNA': 3, 'PosAA': 86, 'PosNA': 89},
                               {'LengthNA': 3, 'PosAA': 87, 'PosNA': 92},
                               {'LengthNA': 3, 'PosAA': 88, 'PosNA': 95},
                               {'LengthNA': 3, 'PosAA': 89, 'PosNA': 98},
                               {'LengthNA': 3, 'PosAA': 90, 'PosNA': 101},
                               {'LengthNA': 3, 'PosAA': 91, 'PosNA': 104},
                               {'LengthNA': 3, 'PosAA': 92, 'PosNA': 107},
                               {'LengthNA': 3, 'PosAA': 93, 'PosNA': 110},
                               {'LengthNA': 3, 'PosAA': 94, 'PosNA': 113},
                               {'LengthNA': 3, 'PosAA': 95, 'PosNA': 116},
                               {'LengthNA': 3, 'PosAA': 96, 'PosNA': 119},
                               {'LengthNA': 3, 'PosAA': 97, 'PosNA': 122},
                               {'LengthNA': 3, 'PosAA': 98, 'PosNA': 125},
                               {'LengthNA': 3, 'PosAA': 99, 'PosNA': 128},
                               {'LengthNA': 3, 'PosAA': 100, 'PosNA': 131},
                               {'LengthNA': 3, 'PosAA': 101, 'PosNA': 134},
                               {'LengthNA': 3, 'PosAA': 102, 'PosNA': 137},
                               {'LengthNA': 3, 'PosAA': 103, 'PosNA': 140},
                               {'LengthNA': 3, 'PosAA': 104, 'PosNA': 143},
                               {'LengthNA': 3, 'PosAA': 105, 'PosNA': 146},
                               {'LengthNA': 3, 'PosAA': 106, 'PosNA': 149},
                               {'LengthNA': 3, 'PosAA': 107, 'PosNA': 152},
                               {'LengthNA': 3, 'PosAA': 108, 'PosNA': 155},
                               {'LengthNA': 3, 'PosAA': 109, 'PosNA': 158},
                               {'LengthNA': 3, 'PosAA': 110, 'PosNA': 161},
                               {'LengthNA': 3, 'PosAA': 111, 'PosNA': 164},
                               {'LengthNA': 3, 'PosAA': 112, 'PosNA': 167},
                               {'LengthNA': 3, 'PosAA': 113, 'PosNA': 170},
                               {'LengthNA': 3, 'PosAA': 114, 'PosNA': 173},
                               {'LengthNA': 3, 'PosAA': 115, 'PosNA': 176},
                               {'LengthNA': 3, 'PosAA': 116, 'PosNA': 179},
                               {'LengthNA': 3, 'PosAA': 117, 'PosNA': 182},
                               {'LengthNA': 3, 'PosAA': 118, 'PosNA': 185},
                               {'LengthNA': 3, 'PosAA': 119, 'PosNA': 188},
                               {'LengthNA': 3, 'PosAA': 120, 'PosNA': 191},
                               {'LengthNA': 3, 'PosAA': 121, 'PosNA': 194},
                               {'LengthNA': 3, 'PosAA': 122, 'PosNA': 197},
                               {'LengthNA': 3, 'PosAA': 123, 'PosNA': 200},
                               {'LengthNA': 3, 'PosAA': 124, 'PosNA': 203},
                               {'LengthNA': 3, 'PosAA': 125, 'PosNA': 206},
                               {'LengthNA': 3, 'PosAA': 126, 'PosNA': 209},
                               {'LengthNA': 3, 'PosAA': 127, 'PosNA': 212},
                               {'LengthNA': 3, 'PosAA': 128, 'PosNA': 215},
                               {'LengthNA': 3, 'PosAA': 129, 'PosNA': 218},
                               {'LengthNA': 3, 'PosAA': 130, 'PosNA': 221},
                               {'LengthNA': 3, 'PosAA': 131, 'PosNA': 224},
                               {'LengthNA': 3, 'PosAA': 132, 'PosNA': 227},
                               {'LengthNA': 3, 'PosAA': 133, 'PosNA': 230},
                               {'LengthNA': 3, 'PosAA': 134, 'PosNA': 233},
                               {'LengthNA': 3, 'PosAA': 135, 'PosNA': 236},
                               {'LengthNA': 3, 'PosAA': 136, 'PosNA': 239},
                               {'LengthNA': 3, 'PosAA': 137, 'PosNA': 242},
                               {'LengthNA': 3, 'PosAA': 138, 'PosNA': 245},
                               {'LengthNA': 3, 'PosAA': 139, 'PosNA': 248},
                               {'LengthNA': 3, 'PosAA': 140, 'PosNA': 251},
                               {'LengthNA': 3, 'PosAA': 141, 'PosNA': 254},
                               {'LengthNA': 3, 'PosAA': 142, 'PosNA': 257},
                               {'LengthNA': 3, 'PosAA': 143, 'PosNA': 260},
                               {'LengthNA': 3, 'PosAA': 144, 'PosNA': 263},
                               {'LengthNA': 3, 'PosAA': 145, 'PosNA': 266},
                               {'LengthNA': 3, 'PosAA': 146, 'PosNA': 269},
                               {'LengthNA': 3, 'PosAA': 147, 'PosNA': 272},
                               {'LengthNA': 3, 'PosAA': 148, 'PosNA': 275},
                               {'LengthNA': 3, 'PosAA': 149, 'PosNA': 278},
                               {'LengthNA': 3, 'PosAA': 150, 'PosNA': 281},
                               {'LengthNA': 3, 'PosAA': 151, 'PosNA': 284},
                               {'LengthNA': 3, 'PosAA': 152, 'PosNA': 287},
                               {'LengthNA': 3, 'PosAA': 153, 'PosNA': 290},
                               {'LengthNA': 3, 'PosAA': 154, 'PosNA': 293},
                               {'LengthNA': 3, 'PosAA': 155, 'PosNA': 296}],
            'FirstAA': 57,
            'FirstNA': 2,
            'Frameshifts': [],
            'LastAA': 155,
            'LastNA': 298,
            'Mutations': [{'AminoAcidText': 'V',
                           'CodonText': 'GTC',
                           'Control': '...',
                           'InsertedAminoAcidsText': '',
                           'InsertedCodonsText': '',
                           'IsDeletion': False,
                           'IsInsertion': False,
                           'IsPartial': False,
                           'NAPosition': 8,
                           'Position': 59,
                           'ReferenceText': 'I'},
                          {'AminoAcidText': 'S',
                           'CodonText': 'AGT',
                           'Control': '...',
                           'InsertedAminoAcidsText': '',
                           'InsertedCodonsText': '',
                           'IsDeletion': False,
                           'IsInsertion': False,
                           'IsPartial': False,
                           'NAPosition': 110,
                           'Position': 93,
                           'ReferenceText': 'N'}],
              'Name': 'plus1',
              'Sequence': 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
             {'AlignedSites': [{'LengthNA': 3, 'PosAA': 56, 'PosNA': 1},
                               {'LengthNA': 3, 'PosAA': 57, 'PosNA': 4},
                               {'LengthNA': 3, 'PosAA': 58, 'PosNA': 7},
                               {'LengthNA': 3, 'PosAA': 59, 'PosNA': 10},
                               {'LengthNA': 3, 'PosAA': 60, 'PosNA': 13},
                               {'LengthNA': 3, 'PosAA': 61, 'PosNA': 16},
                               {'LengthNA': 3, 'PosAA': 62, 'PosNA': 19},
                               {'LengthNA': 3, 'PosAA': 63, 'PosNA': 22},
                               {'LengthNA': 3, 'PosAA': 64, 'PosNA': 25},
                               {'LengthNA': 3, 'PosAA': 65, 'PosNA': 28},
                               {'LengthNA': 3, 'PosAA': 66, 'PosNA': 31},
                               {'LengthNA': 3, 'PosAA': 67, 'PosNA': 34},
                               {'LengthNA': 3, 'PosAA': 68, 'PosNA': 37},
                               {'LengthNA': 3, 'PosAA': 69, 'PosNA': 40},
                               {'LengthNA': 3, 'PosAA': 70, 'PosNA': 43},
                               {'LengthNA': 3, 'PosAA': 71, 'PosNA': 46},
                               {'LengthNA': 3, 'PosAA': 72, 'PosNA': 49},
                               {'LengthNA': 3, 'PosAA': 73, 'PosNA': 52},
                               {'LengthNA': 3, 'PosAA': 74, 'PosNA': 55},
                               {'LengthNA': 3, 'PosAA': 75, 'PosNA': 58},
                               {'LengthNA': 3, 'PosAA': 76, 'PosNA': 61},
                               {'LengthNA': 3, 'PosAA': 77, 'PosNA': 64},
                               {'LengthNA': 3, 'PosAA': 78, 'PosNA': 67},
                               {'LengthNA': 3, 'PosAA': 79, 'PosNA': 70},
                               {'LengthNA': 3, 'PosAA': 80, 'PosNA': 73},
                               {'LengthNA': 3, 'PosAA': 81, 'PosNA': 76},
                               {'LengthNA': 3, 'PosAA': 82, 'PosNA': 79},
                               {'LengthNA': 3, 'PosAA': 83, 'PosNA': 82},
                               {'LengthNA': 3, 'PosAA': 84, 'PosNA': 85},
                               {'LengthNA': 3, 'PosAA': 85, 'PosNA': 88},
                               {'LengthNA': 3, 'PosAA': 86, 'PosNA': 91},
                               {'LengthNA': 3, 'PosAA': 87, 'PosNA': 94},
                               {'LengthNA': 3, 'PosAA': 88, 'PosNA': 97},
                               {'LengthNA': 3, 'PosAA': 89, 'PosNA': 100},
                               {'LengthNA': 3, 'PosAA': 90, 'PosNA': 103},
                               {'LengthNA': 3, 'PosAA': 91, 'PosNA': 106},
                               {'LengthNA': 3, 'PosAA': 92, 'PosNA': 109},
                               {'LengthNA': 3, 'PosAA': 93, 'PosNA': 112},
                               {'LengthNA': 3, 'PosAA': 94, 'PosNA': 115},
                               {'LengthNA': 3, 'PosAA': 95, 'PosNA': 118},
                               {'LengthNA': 3, 'PosAA': 96, 'PosNA': 121},
                               {'LengthNA': 3, 'PosAA': 97, 'PosNA': 124},
                               {'LengthNA': 3, 'PosAA': 98, 'PosNA': 127},
                               {'LengthNA': 3, 'PosAA': 99, 'PosNA': 130},
                               {'LengthNA': 3, 'PosAA': 100, 'PosNA': 133},
                               {'LengthNA': 3, 'PosAA': 101, 'PosNA': 136},
                               {'LengthNA': 3, 'PosAA': 102, 'PosNA': 139},
                               {'LengthNA': 3, 'PosAA': 103, 'PosNA': 142},
                               {'LengthNA': 3, 'PosAA': 104, 'PosNA': 145},
                               {'LengthNA': 3, 'PosAA': 105, 'PosNA': 148},
                               {'LengthNA': 3, 'PosAA': 106, 'PosNA': 151},
                               {'LengthNA': 3, 'PosAA': 107, 'PosNA': 154},
                               {'LengthNA': 3, 'PosAA': 108, 'PosNA': 157},
                               {'LengthNA': 3, 'PosAA': 109, 'PosNA': 160},
                               {'LengthNA': 3, 'PosAA': 110, 'PosNA': 163},
                               {'LengthNA': 3, 'PosAA': 111, 'PosNA': 166},
                               {'LengthNA': 3, 'PosAA': 112, 'PosNA': 169},
                               {'LengthNA': 3, 'PosAA': 113, 'PosNA': 172},
                               {'LengthNA': 3, 'PosAA': 114, 'PosNA': 175},
                               {'LengthNA': 3, 'PosAA': 115, 'PosNA': 178},
                               {'LengthNA': 3, 'PosAA': 116, 'PosNA': 181},
                               {'LengthNA': 3, 'PosAA': 117, 'PosNA': 184},
                               {'LengthNA': 3, 'PosAA': 118, 'PosNA': 187},
                               {'LengthNA': 3, 'PosAA': 119, 'PosNA': 190},
                               {'LengthNA': 3, 'PosAA': 120, 'PosNA': 193},
                               {'LengthNA': 3, 'PosAA': 121, 'PosNA': 196},
                               {'LengthNA': 3, 'PosAA': 122, 'PosNA': 199},
                               {'LengthNA': 3, 'PosAA': 123, 'PosNA': 202},
                               {'LengthNA': 3, 'PosAA': 124, 'PosNA': 205},
                               {'LengthNA': 3, 'PosAA': 125, 'PosNA': 208},
                               {'LengthNA': 3, 'PosAA': 126, 'PosNA': 211},
                               {'LengthNA': 3, 'PosAA': 127, 'PosNA': 214},
                               {'LengthNA': 3, 'PosAA': 128, 'PosNA': 217},
                               {'LengthNA': 3, 'PosAA': 129, 'PosNA': 220},
                               {'LengthNA': 3, 'PosAA': 130, 'PosNA': 223},
                               {'LengthNA': 3, 'PosAA': 131, 'PosNA': 226},
                               {'LengthNA': 3, 'PosAA': 132, 'PosNA': 229},
                               {'LengthNA': 3, 'PosAA': 133, 'PosNA': 232},
                               {'LengthNA': 3, 'PosAA': 134, 'PosNA': 235},
                               {'LengthNA': 3, 'PosAA': 135, 'PosNA': 238},
                               {'LengthNA': 3, 'PosAA': 136, 'PosNA': 241},
                               {'LengthNA': 3, 'PosAA': 137, 'PosNA': 244},
                               {'LengthNA': 3, 'PosAA': 138, 'PosNA': 247},
                               {'LengthNA': 3, 'PosAA': 139, 'PosNA': 250},
                               {'LengthNA': 3, 'PosAA': 140, 'PosNA': 253},
                               {'LengthNA': 3, 'PosAA': 141, 'PosNA': 256},
                               {'LengthNA': 3, 'PosAA': 142, 'PosNA': 259},
                               {'LengthNA': 3, 'PosAA': 143, 'PosNA': 262},
                               {'LengthNA': 3, 'PosAA': 144, 'PosNA': 265},
                               {'LengthNA': 3, 'PosAA': 145, 'PosNA': 268},
                               {'LengthNA': 3, 'PosAA': 146, 'PosNA': 271},
                               {'LengthNA': 3, 'PosAA': 147, 'PosNA': 274},
                               {'LengthNA': 3, 'PosAA': 148, 'PosNA': 277},
                               {'LengthNA': 3, 'PosAA': 149, 'PosNA': 280},
                               {'LengthNA': 3, 'PosAA': 150, 'PosNA': 283},
                               {'LengthNA': 3, 'PosAA': 151, 'PosNA': 286},
                               {'LengthNA': 3, 'PosAA': 152, 'PosNA': 289},
                               {'LengthNA': 3, 'PosAA': 153, 'PosNA': 292},
                               {'LengthNA': 3, 'PosAA': 154, 'PosNA': 295},
                               {'LengthNA': 3, 'PosAA': 155, 'PosNA': 298}],
            'FirstAA': 56,
            'FirstNA': 1,
            'Frameshifts': [],
            'LastAA': 155,
            'LastNA': 300,
            'Mutations': [{'AminoAcidText': 'V',
                            'CodonText': 'GTC',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 10,
                            'Position': 59,
                            'ReferenceText': 'I'},
                            {'AminoAcidText': 'S',
                            'CodonText': 'AGT',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 112,
                            'Position': 93,
                            'ReferenceText': 'N'}],
            'Name': 'plus_codon',
            'Sequence': 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
            {'AlignedSites': [{'LengthNA': 3, 'PosAA': 57, 'PosNA': 1},
                            {'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
                            {'LengthNA': 3, 'PosAA': 59, 'PosNA': 7},
                            {'LengthNA': 2, 'PosAA': 60, 'PosNA': 10},
                            {'LengthNA': 3, 'PosAA': 61, 'PosNA': 12},
                            {'LengthNA': 3, 'PosAA': 62, 'PosNA': 15},
                            {'LengthNA': 3, 'PosAA': 63, 'PosNA': 18},
                            {'LengthNA': 3, 'PosAA': 64, 'PosNA': 21},
                            {'LengthNA': 3, 'PosAA': 65, 'PosNA': 24},
                            {'LengthNA': 3, 'PosAA': 66, 'PosNA': 27},
                            {'LengthNA': 3, 'PosAA': 67, 'PosNA': 30},
                            {'LengthNA': 3, 'PosAA': 68, 'PosNA': 33},
                            {'LengthNA': 3, 'PosAA': 69, 'PosNA': 36},
                            {'LengthNA': 3, 'PosAA': 70, 'PosNA': 39},
                            {'LengthNA': 3, 'PosAA': 71, 'PosNA': 42},
                            {'LengthNA': 3, 'PosAA': 72, 'PosNA': 45},
                            {'LengthNA': 3, 'PosAA': 73, 'PosNA': 48},
                            {'LengthNA': 3, 'PosAA': 74, 'PosNA': 51},
                            {'LengthNA': 3, 'PosAA': 75, 'PosNA': 54},
                            {'LengthNA': 3, 'PosAA': 76, 'PosNA': 57},
                            {'LengthNA': 3, 'PosAA': 77, 'PosNA': 60},
                            {'LengthNA': 3, 'PosAA': 78, 'PosNA': 63},
                            {'LengthNA': 3, 'PosAA': 79, 'PosNA': 66},
                            {'LengthNA': 3, 'PosAA': 80, 'PosNA': 69},
                            {'LengthNA': 3, 'PosAA': 81, 'PosNA': 72},
                            {'LengthNA': 3, 'PosAA': 82, 'PosNA': 75},
                            {'LengthNA': 3, 'PosAA': 83, 'PosNA': 78},
                            {'LengthNA': 3, 'PosAA': 84, 'PosNA': 81},
                            {'LengthNA': 3, 'PosAA': 85, 'PosNA': 84},
                            {'LengthNA': 3, 'PosAA': 86, 'PosNA': 87},
                            {'LengthNA': 3, 'PosAA': 87, 'PosNA': 90},
                            {'LengthNA': 3, 'PosAA': 88, 'PosNA': 93},
                            {'LengthNA': 3, 'PosAA': 89, 'PosNA': 96},
                            {'LengthNA': 3, 'PosAA': 90, 'PosNA': 99},
                            {'LengthNA': 3, 'PosAA': 91, 'PosNA': 102},
                            {'LengthNA': 3, 'PosAA': 92, 'PosNA': 105},
                            {'LengthNA': 3, 'PosAA': 93, 'PosNA': 108},
                            {'LengthNA': 3, 'PosAA': 94, 'PosNA': 111},
                            {'LengthNA': 3, 'PosAA': 95, 'PosNA': 114},
                            {'LengthNA': 3, 'PosAA': 96, 'PosNA': 117},
                            {'LengthNA': 3, 'PosAA': 97, 'PosNA': 120},
                            {'LengthNA': 3, 'PosAA': 98, 'PosNA': 123},
                            {'LengthNA': 3, 'PosAA': 99, 'PosNA': 126},
                            {'LengthNA': 3, 'PosAA': 100, 'PosNA': 129},
                            {'LengthNA': 3, 'PosAA': 101, 'PosNA': 132},
                            {'LengthNA': 3, 'PosAA': 102, 'PosNA': 135},
                            {'LengthNA': 3, 'PosAA': 103, 'PosNA': 138},
                            {'LengthNA': 3, 'PosAA': 104, 'PosNA': 141},
                            {'LengthNA': 3, 'PosAA': 105, 'PosNA': 144},
                            {'LengthNA': 3, 'PosAA': 106, 'PosNA': 147},
                            {'LengthNA': 3, 'PosAA': 107, 'PosNA': 150},
                            {'LengthNA': 3, 'PosAA': 108, 'PosNA': 153},
                            {'LengthNA': 3, 'PosAA': 109, 'PosNA': 156},
                            {'LengthNA': 3, 'PosAA': 110, 'PosNA': 159},
                            {'LengthNA': 3, 'PosAA': 111, 'PosNA': 162},
                            {'LengthNA': 3, 'PosAA': 112, 'PosNA': 165},
                            {'LengthNA': 3, 'PosAA': 113, 'PosNA': 168},
                            {'LengthNA': 3, 'PosAA': 114, 'PosNA': 171},
                            {'LengthNA': 3, 'PosAA': 115, 'PosNA': 174},
                            {'LengthNA': 3, 'PosAA': 116, 'PosNA': 177},
                            {'LengthNA': 3, 'PosAA': 117, 'PosNA': 180},
                            {'LengthNA': 3, 'PosAA': 118, 'PosNA': 183},
                            {'LengthNA': 3, 'PosAA': 119, 'PosNA': 186},
                            {'LengthNA': 3, 'PosAA': 120, 'PosNA': 189},
                            {'LengthNA': 3, 'PosAA': 121, 'PosNA': 192},
                            {'LengthNA': 3, 'PosAA': 122, 'PosNA': 195},
                            {'LengthNA': 3, 'PosAA': 123, 'PosNA': 198},
                            {'LengthNA': 3, 'PosAA': 124, 'PosNA': 201},
                            {'LengthNA': 3, 'PosAA': 125, 'PosNA': 204},
                            {'LengthNA': 3, 'PosAA': 126, 'PosNA': 207},
                            {'LengthNA': 3, 'PosAA': 127, 'PosNA': 210},
                            {'LengthNA': 3, 'PosAA': 128, 'PosNA': 213},
                            {'LengthNA': 3, 'PosAA': 129, 'PosNA': 216},
                            {'LengthNA': 3, 'PosAA': 130, 'PosNA': 219},
                            {'LengthNA': 3, 'PosAA': 131, 'PosNA': 222},
                            {'LengthNA': 3, 'PosAA': 132, 'PosNA': 225},
                            {'LengthNA': 3, 'PosAA': 133, 'PosNA': 228},
                            {'LengthNA': 3, 'PosAA': 134, 'PosNA': 231},
                            {'LengthNA': 3, 'PosAA': 135, 'PosNA': 234},
                            {'LengthNA': 3, 'PosAA': 136, 'PosNA': 237},
                            {'LengthNA': 3, 'PosAA': 137, 'PosNA': 240},
                            {'LengthNA': 3, 'PosAA': 138, 'PosNA': 243},
                            {'LengthNA': 3, 'PosAA': 139, 'PosNA': 246},
                            {'LengthNA': 3, 'PosAA': 140, 'PosNA': 249},
                            {'LengthNA': 3, 'PosAA': 141, 'PosNA': 252},
                            {'LengthNA': 3, 'PosAA': 142, 'PosNA': 255},
                            {'LengthNA': 3, 'PosAA': 143, 'PosNA': 258},
                            {'LengthNA': 3, 'PosAA': 144, 'PosNA': 261},
                            {'LengthNA': 3, 'PosAA': 145, 'PosNA': 264},
                            {'LengthNA': 3, 'PosAA': 146, 'PosNA': 267},
                            {'LengthNA': 3, 'PosAA': 147, 'PosNA': 270},
                            {'LengthNA': 3, 'PosAA': 148, 'PosNA': 273},
                            {'LengthNA': 3, 'PosAA': 149, 'PosNA': 276},
                            {'LengthNA': 3, 'PosAA': 150, 'PosNA': 279},
                            {'LengthNA': 3, 'PosAA': 151, 'PosNA': 282},
                            {'LengthNA': 3, 'PosAA': 152, 'PosNA': 285},
                            {'LengthNA': 3, 'PosAA': 153, 'PosNA': 288},
                            {'LengthNA': 3, 'PosAA': 154, 'PosNA': 291},
                            {'LengthNA': 3, 'PosAA': 155, 'PosNA': 294}],
            'FirstAA': 57,
            'FirstNA': 1,
            'Frameshifts': [{'GapLength': 1,
                            'IsDeletion': True,
                            'IsInsertion': False,
                            'NAPosition': 12,
                            'NucleicAcidsText': '',
                            'Position': 60}],
            'LastAA': 155,
            'LastNA': 296,
            'Mutations': [{'AminoAcidText': 'V',
                            'CodonText': 'GTC',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 7,
                            'Position': 59,
                            'ReferenceText': 'I'},
                            {'AminoAcidText': 'TPAS',
                            'CodonText': ' CT',
                            'Control': '-..',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': True,
                            'NAPosition': 10,
                            'Position': 60,
                            'ReferenceText': 'T'},
                            {'AminoAcidText': 'S',
                            'CodonText': 'AGT',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 108,
                            'Position': 93,
                            'ReferenceText': 'N'}],
            'Name': 'del1_after_3codons',
            'Sequence': 'CCTCAGGTC-CTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'},
            {'AlignedSites': [{'LengthNA': 3, 'PosAA': 57, 'PosNA': 1},
                            {'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
                            {'LengthNA': 6, 'PosAA': 59, 'PosNA': 7},
                            {'LengthNA': 3, 'PosAA': 60, 'PosNA': 13},
                            {'LengthNA': 3, 'PosAA': 61, 'PosNA': 16},
                            {'LengthNA': 3, 'PosAA': 62, 'PosNA': 19},
                            {'LengthNA': 3, 'PosAA': 63, 'PosNA': 22},
                            {'LengthNA': 3, 'PosAA': 64, 'PosNA': 25},
                            {'LengthNA': 3, 'PosAA': 65, 'PosNA': 28},
                            {'LengthNA': 3, 'PosAA': 66, 'PosNA': 31},
                            {'LengthNA': 3, 'PosAA': 67, 'PosNA': 34},
                            {'LengthNA': 3, 'PosAA': 68, 'PosNA': 37},
                            {'LengthNA': 3, 'PosAA': 69, 'PosNA': 40},
                            {'LengthNA': 3, 'PosAA': 70, 'PosNA': 43},
                            {'LengthNA': 3, 'PosAA': 71, 'PosNA': 46},
                            {'LengthNA': 3, 'PosAA': 72, 'PosNA': 49},
                            {'LengthNA': 3, 'PosAA': 73, 'PosNA': 52},
                            {'LengthNA': 3, 'PosAA': 74, 'PosNA': 55},
                            {'LengthNA': 3, 'PosAA': 75, 'PosNA': 58},
                            {'LengthNA': 3, 'PosAA': 76, 'PosNA': 61},
                            {'LengthNA': 3, 'PosAA': 77, 'PosNA': 64},
                            {'LengthNA': 3, 'PosAA': 78, 'PosNA': 67},
                            {'LengthNA': 3, 'PosAA': 79, 'PosNA': 70},
                            {'LengthNA': 3, 'PosAA': 80, 'PosNA': 73},
                            {'LengthNA': 3, 'PosAA': 81, 'PosNA': 76},
                            {'LengthNA': 3, 'PosAA': 82, 'PosNA': 79},
                            {'LengthNA': 3, 'PosAA': 83, 'PosNA': 82},
                            {'LengthNA': 3, 'PosAA': 84, 'PosNA': 85},
                            {'LengthNA': 3, 'PosAA': 85, 'PosNA': 88},
                            {'LengthNA': 3, 'PosAA': 86, 'PosNA': 91},
                            {'LengthNA': 3, 'PosAA': 87, 'PosNA': 94},
                            {'LengthNA': 3, 'PosAA': 88, 'PosNA': 97},
                            {'LengthNA': 3, 'PosAA': 89, 'PosNA': 100},
                            {'LengthNA': 3, 'PosAA': 90, 'PosNA': 103},
                            {'LengthNA': 3, 'PosAA': 91, 'PosNA': 106},
                            {'LengthNA': 3, 'PosAA': 92, 'PosNA': 109},
                            {'LengthNA': 3, 'PosAA': 93, 'PosNA': 112},
                            {'LengthNA': 3, 'PosAA': 94, 'PosNA': 115},
                            {'LengthNA': 3, 'PosAA': 95, 'PosNA': 118},
                            {'LengthNA': 3, 'PosAA': 96, 'PosNA': 121},
                            {'LengthNA': 3, 'PosAA': 97, 'PosNA': 124},
                            {'LengthNA': 3, 'PosAA': 98, 'PosNA': 127},
                            {'LengthNA': 3, 'PosAA': 99, 'PosNA': 130},
                            {'LengthNA': 3, 'PosAA': 100, 'PosNA': 133},
                            {'LengthNA': 3, 'PosAA': 101, 'PosNA': 136},
                            {'LengthNA': 3, 'PosAA': 102, 'PosNA': 139},
                            {'LengthNA': 3, 'PosAA': 103, 'PosNA': 142},
                            {'LengthNA': 3, 'PosAA': 104, 'PosNA': 145},
                            {'LengthNA': 3, 'PosAA': 105, 'PosNA': 148},
                            {'LengthNA': 3, 'PosAA': 106, 'PosNA': 151},
                            {'LengthNA': 3, 'PosAA': 107, 'PosNA': 154},
                            {'LengthNA': 3, 'PosAA': 108, 'PosNA': 157},
                            {'LengthNA': 3, 'PosAA': 109, 'PosNA': 160},
                            {'LengthNA': 3, 'PosAA': 110, 'PosNA': 163},
                            {'LengthNA': 3, 'PosAA': 111, 'PosNA': 166},
                            {'LengthNA': 3, 'PosAA': 112, 'PosNA': 169},
                            {'LengthNA': 3, 'PosAA': 113, 'PosNA': 172},
                            {'LengthNA': 3, 'PosAA': 114, 'PosNA': 175},
                            {'LengthNA': 3, 'PosAA': 115, 'PosNA': 178},
                            {'LengthNA': 3, 'PosAA': 116, 'PosNA': 181},
                            {'LengthNA': 3, 'PosAA': 117, 'PosNA': 184},
                            {'LengthNA': 3, 'PosAA': 118, 'PosNA': 187},
                            {'LengthNA': 3, 'PosAA': 119, 'PosNA': 190},
                            {'LengthNA': 3, 'PosAA': 120, 'PosNA': 193},
                            {'LengthNA': 3, 'PosAA': 121, 'PosNA': 196},
                            {'LengthNA': 3, 'PosAA': 122, 'PosNA': 199},
                            {'LengthNA': 3, 'PosAA': 123, 'PosNA': 202},
                            {'LengthNA': 3, 'PosAA': 124, 'PosNA': 205},
                            {'LengthNA': 3, 'PosAA': 125, 'PosNA': 208},
                            {'LengthNA': 3, 'PosAA': 126, 'PosNA': 211},
                            {'LengthNA': 3, 'PosAA': 127, 'PosNA': 214},
                            {'LengthNA': 3, 'PosAA': 128, 'PosNA': 217},
                            {'LengthNA': 3, 'PosAA': 129, 'PosNA': 220},
                            {'LengthNA': 3, 'PosAA': 130, 'PosNA': 223},
                            {'LengthNA': 3, 'PosAA': 131, 'PosNA': 226},
                            {'LengthNA': 3, 'PosAA': 132, 'PosNA': 229},
                            {'LengthNA': 3, 'PosAA': 133, 'PosNA': 232},
                            {'LengthNA': 3, 'PosAA': 134, 'PosNA': 235},
                            {'LengthNA': 3, 'PosAA': 135, 'PosNA': 238},
                            {'LengthNA': 3, 'PosAA': 136, 'PosNA': 241},
                            {'LengthNA': 3, 'PosAA': 137, 'PosNA': 244},
                            {'LengthNA': 3, 'PosAA': 138, 'PosNA': 247},
                            {'LengthNA': 3, 'PosAA': 139, 'PosNA': 250},
                            {'LengthNA': 3, 'PosAA': 140, 'PosNA': 253},
                            {'LengthNA': 3, 'PosAA': 141, 'PosNA': 256},
                            {'LengthNA': 3, 'PosAA': 142, 'PosNA': 259},
                            {'LengthNA': 3, 'PosAA': 143, 'PosNA': 262},
                            {'LengthNA': 3, 'PosAA': 144, 'PosNA': 265},
                            {'LengthNA': 3, 'PosAA': 145, 'PosNA': 268},
                            {'LengthNA': 3, 'PosAA': 146, 'PosNA': 271},
                            {'LengthNA': 3, 'PosAA': 147, 'PosNA': 274},
                            {'LengthNA': 3, 'PosAA': 148, 'PosNA': 277},
                            {'LengthNA': 3, 'PosAA': 149, 'PosNA': 280},
                            {'LengthNA': 3, 'PosAA': 150, 'PosNA': 283},
                            {'LengthNA': 3, 'PosAA': 151, 'PosNA': 286},
                            {'LengthNA': 3, 'PosAA': 152, 'PosNA': 289},
                            {'LengthNA': 3, 'PosAA': 153, 'PosNA': 292},
                            {'LengthNA': 3, 'PosAA': 154, 'PosNA': 295},
                            {'LengthNA': 3, 'PosAA': 155, 'PosNA': 298}],
            'FirstAA': 57,
            'FirstNA': 1,
            'Frameshifts': [],
            'LastAA': 155,
            'LastNA': 300,
            'Mutations': [{'AminoAcidText': 'V',
                            'CodonText': 'GTC',
                            'Control': '...+++',
                            'InsertedAminoAcidsText': 'K',
                            'InsertedCodonsText': 'AAA',
                            'IsDeletion': False,
                            'IsInsertion': True,
                            'IsPartial': False,
                            'NAPosition': 7,
                            'Position': 59,
                            'ReferenceText': 'I'},
                            {'AminoAcidText': 'S',
                            'CodonText': 'AGT',
                            'Control': '...',
                            'InsertedAminoAcidsText': '',
                            'InsertedCodonsText': '',
                            'IsDeletion': False,
                            'IsInsertion': False,
                            'IsPartial': False,
                            'NAPosition': 112,
                            'Position': 93,
                            'ReferenceText': 'N'}],
            'Name': 'insAAA_after_3codons',
            'Sequence': 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAAT'}]

        res_records = self.aligner.align_file(filename, program='nuc')

        self.assertEqual(exp_records, res_records)

    def testAlignFilePost(self):
        # Setting params
        filename = r'tests/hxb2-pr.fa'

        exp_records = \
            [
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 3,
                            "PosAA": 57,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 4
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 7
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 10
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 13
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 16
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 19
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 22
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 25
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 28
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 31
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 34
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 37
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 40
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 43
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 46
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 49
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 52
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 55
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 58
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 61
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 64
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 67
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 70
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 73
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 76
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 79
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 82
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 85
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 88
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 91
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 94
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 97
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 100
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 103
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 106
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 109
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 112
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 115
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 118
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 121
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 124
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 127
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 130
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 133
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 136
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 139
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 142
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 145
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 148
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 151
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 154
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 157
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 160
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 163
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 166
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 169
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 172
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 175
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 178
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 181
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 184
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 187
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 190
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 193
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 196
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 199
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 202
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 205
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 208
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 211
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 214
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 217
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 220
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 223
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 226
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 229
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 232
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 235
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 238
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 241
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 244
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 247
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 250
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 253
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 256
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 259
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 262
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 265
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 268
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 271
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 274
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 277
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 280
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 283
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 286
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 289
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 292
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 295
                        }
                    ],
                    "FirstAA": 57,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 297,
                    "Mutations": [
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "HXB2-PR",
                    "Sequence": "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 2,
                            "PosAA": 57,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 3
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 6
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 9
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 12
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 15
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 18
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 21
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 24
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 27
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 30
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 33
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 36
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 39
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 42
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 45
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 48
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 51
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 54
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 57
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 60
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 63
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 66
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 69
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 72
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 75
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 78
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 81
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 84
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 87
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 90
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 93
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 96
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 99
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 102
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 105
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 108
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 111
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 114
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 117
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 120
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 123
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 126
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 129
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 132
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 135
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 138
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 141
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 144
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 147
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 150
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 153
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 156
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 159
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 162
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 165
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 168
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 171
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 174
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 177
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 180
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 183
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 186
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 189
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 192
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 195
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 198
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 201
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 204
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 207
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 210
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 213
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 216
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 219
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 222
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 225
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 228
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 231
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 234
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 237
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 240
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 243
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 246
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 249
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 252
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 255
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 258
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 261
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 264
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 267
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 270
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 273
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 276
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 279
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 282
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 285
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 288
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 291
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 294
                        }
                    ],
                    "FirstAA": 57,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 296,
                    "Mutations": [
                        {
                            "AminoAcidText": "X",
                            "CodonText": "CT-",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 57,
                            "RefCodonText": "CCT",
                            "ReferenceText": "P"
                        },
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "shift1",
                    "Sequence": "CTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT\n"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 1,
                            "PosAA": 57,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 2
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 5
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 8
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 11
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 14
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 17
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 20
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 23
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 26
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 29
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 32
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 35
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 38
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 41
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 44
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 47
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 50
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 53
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 56
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 59
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 62
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 65
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 68
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 71
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 74
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 77
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 80
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 83
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 86
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 89
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 92
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 95
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 98
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 101
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 104
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 107
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 110
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 113
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 116
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 119
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 122
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 125
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 128
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 131
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 134
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 137
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 140
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 143
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 146
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 149
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 152
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 155
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 158
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 161
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 164
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 167
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 170
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 173
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 176
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 179
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 182
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 185
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 188
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 191
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 194
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 197
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 200
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 203
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 206
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 209
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 212
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 215
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 218
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 221
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 224
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 227
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 230
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 233
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 236
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 239
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 242
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 245
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 248
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 251
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 254
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 257
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 260
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 263
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 266
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 269
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 272
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 275
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 278
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 281
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 284
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 287
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 290
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 293
                        }
                    ],
                    "FirstAA": 57,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 295,
                    "Mutations": [
                        {
                            "AminoAcidText": "X",
                            "CodonText": "T--",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 57,
                            "RefCodonText": "CCT",
                            "ReferenceText": "P"
                        },
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "shift2",
                    "Sequence": "TCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT\n"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 1,
                            "PosAA": 56,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 57,
                            "PosNA": 2
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 5
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 8
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 11
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 14
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 17
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 20
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 23
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 26
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 29
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 32
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 35
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 38
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 41
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 44
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 47
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 50
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 53
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 56
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 59
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 62
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 65
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 68
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 71
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 74
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 77
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 80
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 83
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 86
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 89
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 92
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 95
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 98
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 101
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 104
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 107
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 110
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 113
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 116
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 119
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 122
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 125
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 128
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 131
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 134
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 137
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 140
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 143
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 146
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 149
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 152
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 155
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 158
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 161
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 164
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 167
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 170
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 173
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 176
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 179
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 182
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 185
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 188
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 191
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 194
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 197
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 200
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 203
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 206
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 209
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 212
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 215
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 218
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 221
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 224
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 227
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 230
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 233
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 236
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 239
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 242
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 245
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 248
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 251
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 254
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 257
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 260
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 263
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 266
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 269
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 272
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 275
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 278
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 281
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 284
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 287
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 290
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 293
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 296
                        }
                    ],
                    "FirstAA": 56,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 298,
                    "Mutations": [
                        {
                            "AminoAcidText": "X",
                            "CodonText": "C--",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 56,
                            "RefCodonText": "TTC",
                            "ReferenceText": "F"
                        },
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "plus1",
                    "Sequence": "TCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT\n"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 3,
                            "PosAA": 56,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 57,
                            "PosNA": 4
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 7
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 10
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 13
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 16
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 19
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 22
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 25
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 28
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 31
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 34
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 37
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 40
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 43
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 46
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 49
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 52
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 55
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 58
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 61
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 64
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 67
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 70
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 73
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 76
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 79
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 82
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 85
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 88
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 91
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 94
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 97
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 100
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 103
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 106
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 109
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 112
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 115
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 118
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 121
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 124
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 127
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 130
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 133
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 136
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 139
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 142
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 145
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 148
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 151
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 154
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 157
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 160
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 163
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 166
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 169
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 172
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 175
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 178
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 181
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 184
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 187
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 190
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 193
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 196
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 199
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 202
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 205
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 208
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 211
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 214
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 217
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 220
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 223
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 226
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 229
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 232
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 235
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 238
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 241
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 244
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 247
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 250
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 253
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 256
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 259
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 262
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 265
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 268
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 271
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 274
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 277
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 280
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 283
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 286
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 289
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 292
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 295
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 298
                        }
                    ],
                    "FirstAA": 56,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 300,
                    "Mutations": [
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "plus_codon",
                    "Sequence": "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 2,
                            "PosAA": 57,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 3
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 6
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 9
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 12
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 15
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 18
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 21
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 24
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 27
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 30
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 33
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 36
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 39
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 42
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 45
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 48
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 51
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 54
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 57
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 60
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 63
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 66
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 69
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 72
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 75
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 78
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 81
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 84
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 87
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 90
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 93
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 96
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 99
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 102
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 105
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 108
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 111
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 114
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 117
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 120
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 123
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 126
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 129
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 132
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 135
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 138
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 141
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 144
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 147
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 150
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 153
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 156
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 159
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 162
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 165
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 168
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 171
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 174
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 177
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 180
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 183
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 186
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 189
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 192
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 195
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 198
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 201
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 204
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 207
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 210
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 213
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 216
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 219
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 222
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 225
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 228
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 231
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 234
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 237
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 240
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 243
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 246
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 249
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 252
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 255
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 258
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 261
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 264
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 267
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 270
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 273
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 276
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 279
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 282
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 285
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 288
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 291
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 294
                        }
                    ],
                    "FirstAA": 57,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 296,
                    "Mutations": [
                        {
                            "AminoAcidText": "X",
                            "CodonText": "CC-",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 57,
                            "RefCodonText": "CCT",
                            "ReferenceText": "P"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "TCA",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 58,
                            "RefCodonText": "CAG",
                            "ReferenceText": "Q"
                        },
                        {
                            "AminoAcidText": "G",
                            "CodonText": "GGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "P",
                            "CodonText": "CCT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 60,
                            "RefCodonText": "ACT",
                            "ReferenceText": "T"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "del1_after_3codons",
                    "Sequence": "CCTCAGGTCCTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT\n"
                },
                {
                    "AlignedSites": [
                        {
                            "LengthNA": 3,
                            "PosAA": 56,
                            "PosNA": 1
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 57,
                            "PosNA": 4
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 58,
                            "PosNA": 7
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 59,
                            "PosNA": 10
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 60,
                            "PosNA": 13
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 61,
                            "PosNA": 16
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 62,
                            "PosNA": 19
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 63,
                            "PosNA": 22
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 64,
                            "PosNA": 25
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 65,
                            "PosNA": 28
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 66,
                            "PosNA": 31
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 67,
                            "PosNA": 34
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 68,
                            "PosNA": 37
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 69,
                            "PosNA": 40
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 70,
                            "PosNA": 43
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 71,
                            "PosNA": 46
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 72,
                            "PosNA": 49
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 73,
                            "PosNA": 52
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 74,
                            "PosNA": 55
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 75,
                            "PosNA": 58
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 76,
                            "PosNA": 61
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 77,
                            "PosNA": 64
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 78,
                            "PosNA": 67
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 79,
                            "PosNA": 70
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 80,
                            "PosNA": 73
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 81,
                            "PosNA": 76
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 82,
                            "PosNA": 79
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 83,
                            "PosNA": 82
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 84,
                            "PosNA": 85
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 85,
                            "PosNA": 88
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 86,
                            "PosNA": 91
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 87,
                            "PosNA": 94
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 88,
                            "PosNA": 97
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 89,
                            "PosNA": 100
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 90,
                            "PosNA": 103
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 91,
                            "PosNA": 106
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 92,
                            "PosNA": 109
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 93,
                            "PosNA": 112
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 94,
                            "PosNA": 115
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 95,
                            "PosNA": 118
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 96,
                            "PosNA": 121
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 97,
                            "PosNA": 124
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 98,
                            "PosNA": 127
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 99,
                            "PosNA": 130
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 100,
                            "PosNA": 133
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 101,
                            "PosNA": 136
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 102,
                            "PosNA": 139
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 103,
                            "PosNA": 142
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 104,
                            "PosNA": 145
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 105,
                            "PosNA": 148
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 106,
                            "PosNA": 151
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 107,
                            "PosNA": 154
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 108,
                            "PosNA": 157
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 109,
                            "PosNA": 160
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 110,
                            "PosNA": 163
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 111,
                            "PosNA": 166
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 112,
                            "PosNA": 169
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 113,
                            "PosNA": 172
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 114,
                            "PosNA": 175
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 115,
                            "PosNA": 178
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 116,
                            "PosNA": 181
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 117,
                            "PosNA": 184
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 118,
                            "PosNA": 187
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 119,
                            "PosNA": 190
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 120,
                            "PosNA": 193
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 121,
                            "PosNA": 196
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 122,
                            "PosNA": 199
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 123,
                            "PosNA": 202
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 124,
                            "PosNA": 205
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 125,
                            "PosNA": 208
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 126,
                            "PosNA": 211
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 127,
                            "PosNA": 214
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 128,
                            "PosNA": 217
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 129,
                            "PosNA": 220
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 130,
                            "PosNA": 223
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 131,
                            "PosNA": 226
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 132,
                            "PosNA": 229
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 133,
                            "PosNA": 232
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 134,
                            "PosNA": 235
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 135,
                            "PosNA": 238
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 136,
                            "PosNA": 241
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 137,
                            "PosNA": 244
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 138,
                            "PosNA": 247
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 139,
                            "PosNA": 250
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 140,
                            "PosNA": 253
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 141,
                            "PosNA": 256
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 142,
                            "PosNA": 259
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 143,
                            "PosNA": 262
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 144,
                            "PosNA": 265
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 145,
                            "PosNA": 268
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 146,
                            "PosNA": 271
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 147,
                            "PosNA": 274
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 148,
                            "PosNA": 277
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 149,
                            "PosNA": 280
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 150,
                            "PosNA": 283
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 151,
                            "PosNA": 286
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 152,
                            "PosNA": 289
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 153,
                            "PosNA": 292
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 154,
                            "PosNA": 295
                        },
                        {
                            "LengthNA": 3,
                            "PosAA": 155,
                            "PosNA": 298
                        }
                    ],
                    "FirstAA": 56,
                    "FirstNA": 1,
                    "FrameShifts": [],
                    "LastAA": 155,
                    "LastNA": 300,
                    "Mutations": [
                        {
                            "AminoAcidText": "P",
                            "CodonText": "CCT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 56,
                            "RefCodonText": "TTC",
                            "ReferenceText": "F"
                        },
                        {
                            "AminoAcidText": "Q",
                            "CodonText": "CAG",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 57,
                            "RefCodonText": "CCT",
                            "ReferenceText": "P"
                        },
                        {
                            "AminoAcidText": "V",
                            "CodonText": "GTC",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 58,
                            "RefCodonText": "CAG",
                            "ReferenceText": "Q"
                        },
                        {
                            "AminoAcidText": "K",
                            "CodonText": "AAA",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 59,
                            "RefCodonText": "ATC",
                            "ReferenceText": "I"
                        },
                        {
                            "AminoAcidText": "S",
                            "CodonText": "AGT",
                            "InsertedCodonsText": "",
                            "IsDeletion": False,
                            "IsInsertion": False,
                            "Position": 93,
                            "RefCodonText": "AAT",
                            "ReferenceText": "N"
                        }
                    ],
                    "Name": "insAAA_after_3codons",
                    "Sequence": ""
                }
            ]

        res_records = self.aligner.align_file(filename)

        self.assertEqual(exp_records, res_records)

    def testCreateGeneMap(self):
        # Setting params
        exp_map = {'IN': (715, 1003),
                   'PR': (56, 154),
                   'RT': (155, 714)}
        res_map = self.aligner.create_gene_map()

        self.assertEqual(exp_map, res_map)        

    def testGetGenes(self):
        # Setting params
        pol_aligned_sites = \
            [{'LengthNA': 3, 'PosAA': 57, 'PosNA': 1},
             {'LengthNA': 3, 'PosAA': 58, 'PosNA': 4},
             {'LengthNA': 3, 'PosAA': 59, 'PosNA': 7},
             {'LengthNA': 3, 'PosAA': 60, 'PosNA': 10},
             {'LengthNA': 3, 'PosAA': 61, 'PosNA': 13},
             {'LengthNA': 3, 'PosAA': 62, 'PosNA': 16},
             {'LengthNA': 3, 'PosAA': 63, 'PosNA': 19},
             {'LengthNA': 3, 'PosAA': 64, 'PosNA': 22},
             {'LengthNA': 3, 'PosAA': 65, 'PosNA': 25},
             {'LengthNA': 3, 'PosAA': 66, 'PosNA': 28},
             {'LengthNA': 3, 'PosAA': 67, 'PosNA': 31},
             {'LengthNA': 3, 'PosAA': 68, 'PosNA': 34},
             {'LengthNA': 3, 'PosAA': 69, 'PosNA': 37},
             {'LengthNA': 3, 'PosAA': 70, 'PosNA': 40},
             {'LengthNA': 3, 'PosAA': 71, 'PosNA': 43},
             {'LengthNA': 3, 'PosAA': 72, 'PosNA': 46},
             {'LengthNA': 3, 'PosAA': 73, 'PosNA': 49},
             {'LengthNA': 3, 'PosAA': 74, 'PosNA': 52},
             {'LengthNA': 3, 'PosAA': 75, 'PosNA': 55},
             {'LengthNA': 3, 'PosAA': 76, 'PosNA': 58},
             {'LengthNA': 3, 'PosAA': 77, 'PosNA': 61},
             {'LengthNA': 3, 'PosAA': 78, 'PosNA': 64},
             {'LengthNA': 3, 'PosAA': 79, 'PosNA': 67},
             {'LengthNA': 3, 'PosAA': 80, 'PosNA': 70},
             {'LengthNA': 3, 'PosAA': 81, 'PosNA': 73},
             {'LengthNA': 3, 'PosAA': 82, 'PosNA': 76},
             {'LengthNA': 3, 'PosAA': 83, 'PosNA': 79},
             {'LengthNA': 3, 'PosAA': 84, 'PosNA': 82},
             {'LengthNA': 3, 'PosAA': 85, 'PosNA': 85},
             {'LengthNA': 3, 'PosAA': 86, 'PosNA': 88},
             {'LengthNA': 3, 'PosAA': 87, 'PosNA': 91},
             {'LengthNA': 3, 'PosAA': 88, 'PosNA': 94},
             {'LengthNA': 3, 'PosAA': 89, 'PosNA': 97},
             {'LengthNA': 3, 'PosAA': 90, 'PosNA': 100},
             {'LengthNA': 3, 'PosAA': 91, 'PosNA': 103},
             {'LengthNA': 3, 'PosAA': 92, 'PosNA': 106},
             {'LengthNA': 3, 'PosAA': 93, 'PosNA': 109},
             {'LengthNA': 3, 'PosAA': 94, 'PosNA': 112},
             {'LengthNA': 3, 'PosAA': 95, 'PosNA': 115},
             {'LengthNA': 3, 'PosAA': 96, 'PosNA': 118},
             {'LengthNA': 3, 'PosAA': 97, 'PosNA': 121},
             {'LengthNA': 3, 'PosAA': 98, 'PosNA': 124},
             {'LengthNA': 3, 'PosAA': 99, 'PosNA': 127},
             {'LengthNA': 3, 'PosAA': 100, 'PosNA': 130},
             {'LengthNA': 3, 'PosAA': 101, 'PosNA': 133},
             {'LengthNA': 3, 'PosAA': 102, 'PosNA': 136},
             {'LengthNA': 3, 'PosAA': 103, 'PosNA': 139},
             {'LengthNA': 3, 'PosAA': 104, 'PosNA': 142},
             {'LengthNA': 3, 'PosAA': 105, 'PosNA': 145},
             {'LengthNA': 3, 'PosAA': 106, 'PosNA': 148},
             {'LengthNA': 3, 'PosAA': 107, 'PosNA': 151},
             {'LengthNA': 3, 'PosAA': 108, 'PosNA': 154},
             {'LengthNA': 3, 'PosAA': 109, 'PosNA': 157},
             {'LengthNA': 3, 'PosAA': 110, 'PosNA': 160},
             {'LengthNA': 3, 'PosAA': 111, 'PosNA': 163},
             {'LengthNA': 3, 'PosAA': 112, 'PosNA': 166},
             {'LengthNA': 3, 'PosAA': 113, 'PosNA': 169},
             {'LengthNA': 3, 'PosAA': 114, 'PosNA': 172},
             {'LengthNA': 3, 'PosAA': 115, 'PosNA': 175},
             {'LengthNA': 3, 'PosAA': 116, 'PosNA': 178},
             {'LengthNA': 3, 'PosAA': 117, 'PosNA': 181},
             {'LengthNA': 3, 'PosAA': 118, 'PosNA': 184},
             {'LengthNA': 3, 'PosAA': 119, 'PosNA': 187},
             {'LengthNA': 3, 'PosAA': 120, 'PosNA': 190},
             {'LengthNA': 3, 'PosAA': 121, 'PosNA': 193},
             {'LengthNA': 3, 'PosAA': 122, 'PosNA': 196},
             {'LengthNA': 3, 'PosAA': 123, 'PosNA': 199},
             {'LengthNA': 3, 'PosAA': 124, 'PosNA': 202},
             {'LengthNA': 3, 'PosAA': 125, 'PosNA': 205},
             {'LengthNA': 3, 'PosAA': 126, 'PosNA': 208},
             {'LengthNA': 3, 'PosAA': 127, 'PosNA': 211},
             {'LengthNA': 3, 'PosAA': 128, 'PosNA': 214},
             {'LengthNA': 3, 'PosAA': 129, 'PosNA': 217},
             {'LengthNA': 3, 'PosAA': 130, 'PosNA': 220},
             {'LengthNA': 3, 'PosAA': 131, 'PosNA': 223},
             {'LengthNA': 3, 'PosAA': 132, 'PosNA': 226},
             {'LengthNA': 3, 'PosAA': 133, 'PosNA': 229},
             {'LengthNA': 3, 'PosAA': 134, 'PosNA': 232},
             {'LengthNA': 3, 'PosAA': 135, 'PosNA': 235},
             {'LengthNA': 3, 'PosAA': 136, 'PosNA': 238},
             {'LengthNA': 3, 'PosAA': 137, 'PosNA': 241},
             {'LengthNA': 3, 'PosAA': 138, 'PosNA': 244},
             {'LengthNA': 3, 'PosAA': 139, 'PosNA': 247},
             {'LengthNA': 3, 'PosAA': 140, 'PosNA': 250},
             {'LengthNA': 3, 'PosAA': 141, 'PosNA': 253},
             {'LengthNA': 3, 'PosAA': 142, 'PosNA': 256},
             {'LengthNA': 3, 'PosAA': 143, 'PosNA': 259},
             {'LengthNA': 3, 'PosAA': 144, 'PosNA': 262},
             {'LengthNA': 3, 'PosAA': 145, 'PosNA': 265},
             {'LengthNA': 3, 'PosAA': 146, 'PosNA': 268},
             {'LengthNA': 3, 'PosAA': 147, 'PosNA': 271},
             {'LengthNA': 3, 'PosAA': 148, 'PosNA': 274},
             {'LengthNA': 3, 'PosAA': 149, 'PosNA': 277},
             {'LengthNA': 3, 'PosAA': 150, 'PosNA': 280},
             {'LengthNA': 3, 'PosAA': 151, 'PosNA': 283},
             {'LengthNA': 3, 'PosAA': 152, 'PosNA': 286},
             {'LengthNA': 3, 'PosAA': 153, 'PosNA': 289},
             {'LengthNA': 3, 'PosAA': 154, 'PosNA': 292},
             {'LengthNA': 3, 'PosAA': 155, 'PosNA': 295}]
        pol_first_aa = 57
        pol_last_aa = 155

        exp_genes = [('PR', 1, 99, 1, 294)]
        res_genes = self.aligner.get_genes(pol_aligned_sites, pol_first_aa, pol_last_aa)

        self.assertEqual(exp_genes, res_genes)

        # Setting params
        pol_first_aa = 123
        pol_last_aa = 38

        exp_genes = []
        res_genes = self.aligner.get_genes(pol_aligned_sites, pol_first_aa, pol_last_aa)

        self.assertEqual(exp_genes, res_genes)

    def testGetMutations(self):
        # Setting params
        self.maxDiff = None
        records = self.aligner.align_file(r'tests/hxb2-pr.fa', program='nuc')
        
        exp_sequence_headers = ['HXB2-PR', 'shift1', 'shift2', 'plus1', 'plus_codon', 'del1_after_3codons', 'insAAA_after_3codons']
        exp_file_genes = [[('PR', 1, 99, 1, 294)], [('PR', 2, 99, 3, 293)], [('PR', 2, 99, 2, 292)], [('PR', 1, 99, 2, 295)], [('PR', 1, 99, 1, 297)], [('PR', 1, 99, 1, 293)], [('PR', 1, 99, 1, 297)]]  
        exp_file_mutations = [[{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 4: ('T', 'TPAS', 'X'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}]]
        exp_file_trims = [[(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
        exp_subtypes = ['', '', '', '', '', '', '']

        res_sequence_headers, res_file_genes, res_file_mutations, \
            res_file_trims, res_subtypes = self.aligner.get_mutations(records)

        self.assertEqual(exp_sequence_headers, res_sequence_headers)
        self.assertEqual(exp_file_genes, res_file_genes)
        self.assertEqual(exp_file_mutations, res_file_mutations)
        self.assertEqual(exp_file_trims, res_file_trims)
        self.assertEqual(exp_subtypes, res_subtypes)

        # Setting params
        exp_sequence_headers = ['HXB2-PR', 'shift1', 'shift2', 'plus1', 'plus_codon', 'del1_after_3codons', 'insAAA_after_3codons']
        exp_file_genes = [[('PR', 1, 99, 1, 294)], [('PR', 2, 99, 3, 293)], [('PR', 2, 99, 2, 292)], [('PR', 1, 99, 2, 295)], [('PR', 1, 99, 1, 297)], [('PR', 1, 99, 1, 293)], [('PR', 1, 99, 1, 297)]]  
        exp_file_mutations = [[{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 4: ('T', 'TPAS', 'X'), 37: ('N', 'S', 'S')}],
                              [{3: ('I', 'V', 'V'), 37: ('N', 'S', 'S')}]]
        exp_file_trims = [[(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)], [(0, 0)]]
        exp_subtypes = ['B', 'B', 'B', 'B', 'B', 'B', 'B']

        res_sequence_headers, res_file_genes, res_file_mutations, \
            res_file_trims, res_subtypes = self.aligner.get_mutations(records, do_subtype=True)
        
        self.assertEqual(exp_sequence_headers, res_sequence_headers)
        self.assertEqual(exp_file_genes, res_file_genes)
        self.assertEqual(exp_file_mutations, res_file_mutations)
        self.assertEqual(exp_file_trims, res_file_trims)
        self.assertEqual(exp_subtypes, res_subtypes)

    def testTranslateNATriplet(self):
        # Setting params
        triplet = 'YTD'

        exp_translation = 'LF'
        res_translation = self.aligner.translate_na_triplet(triplet)

        self.assertEqual(exp_translation, res_translation)

        # Setting params
        triplet = ''    # len(triplet) == 0

        exp_translation = '-'
        res_translation = self.aligner.translate_na_triplet(triplet)

        self.assertEqual(exp_translation, res_translation)

        # Setting params
        triplet = '~TD' # '~' in triplet

        exp_translation = 'X'
        res_translation = self.aligner.translate_na_triplet(triplet)

        self.assertEqual(exp_translation, res_translation)

        # Setting params
        triplet = 'YTDF' # len(triplet) != 3

        exp_translation = 'X'
        res_translation = self.aligner.translate_na_triplet(triplet)

        self.assertEqual(exp_translation, res_translation)

    def testGenerateTable(self):
        # Setting params
        exp_triplet_table = \
            {'AAA': 'K',
             'AAB': 'NK',
             'AAC': 'N',
             'AAD': 'KN',
             'AAG': 'K',
             'AAH': 'KN',
             'AAK': 'KN',
             'AAM': 'KN',
             'AAN': 'KN',
             'AAR': 'K',
             'AAS': 'NK',
             'AAT': 'N',
             'AAV': 'KN',
             'AAW': 'KN',
             'AAY': 'N',
             'ABA': 'TRI',
             'ABB': 'X',
             'ABC': 'TSI',
             'ABD': 'X',
             'ABG': 'TRM',
             'ABH': 'TRSI',
             'ABK': 'X',
             'ABM': 'TRSI',
             'ABN': 'X',
             'ABR': 'TRIM',
             'ABS': 'X',
             'ABT': 'TSI',
             'ABV': 'X',
             'ABW': 'TRSI',
             'ABY': 'TSI',
             'ACA': 'T',
             'ACB': 'T',
             'ACC': 'T',
             'ACD': 'T',
             'ACG': 'T',
             'ACH': 'T',
             'ACK': 'T',
             'ACM': 'T',
             'ACN': 'T',
             'ACR': 'T',
             'ACS': 'T',
             'ACT': 'T',
             'ACV': 'T',
             'ACW': 'T',
             'ACY': 'T',
             'ADA': 'KRI',
             'ADB': 'X',
             'ADC': 'NSI',
             'ADD': 'X',
             'ADG': 'KRM',
             'ADH': 'X',
             'ADK': 'X',
             'ADM': 'X',
             'ADN': 'X',
             'ADR': 'KRIM',
             'ADS': 'X',
             'ADT': 'NSI',
             'ADV': 'X',
             'ADW': 'X',
             'ADY': 'NSI',
             'AGA': 'R',
             'AGB': 'SR',
             'AGC': 'S',
             'AGD': 'RS',
             'AGG': 'R',
             'AGH': 'RS',
             'AGK': 'RS',
             'AGM': 'RS',
             'AGN': 'RS',
             'AGR': 'R',
             'AGS': 'SR',
             'AGT': 'S',
             'AGV': 'RS',
             'AGW': 'RS',
             'AGY': 'S',
             'AHA': 'KTI',
             'AHB': 'X',
             'AHC': 'NTI',
             'AHD': 'X',
             'AHG': 'KTM',
             'AHH': 'KNTI',
             'AHK': 'X',
             'AHM': 'KNTI',
             'AHN': 'X',
             'AHR': 'KTIM',
             'AHS': 'X',
             'AHT': 'NTI',
             'AHV': 'X',
             'AHW': 'KNTI',
             'AHY': 'NTI',
             'AKA': 'RI',
             'AKB': 'SRIM',
             'AKC': 'SI',
             'AKD': 'RSIM',
             'AKG': 'RM',
             'AKH': 'RSI',
             'AKK': 'RSMI',
             'AKM': 'RSI',
             'AKN': 'RSIM',
             'AKR': 'RIM',
             'AKS': 'SRIM',
             'AKT': 'SI',
             'AKV': 'RSIM',
             'AKW': 'RSI',
             'AKY': 'SI',
             'AMA': 'KT',
             'AMB': 'NKT',
             'AMC': 'NT',
             'AMD': 'KNT',
             'AMG': 'KT',
             'AMH': 'KNT',
             'AMK': 'KNT',
             'AMM': 'KNT',
             'AMN': 'KNT',
             'AMR': 'KT',
             'AMS': 'NKT',
             'AMT': 'NT',
             'AMV': 'KNT',
             'AMW': 'KNT',
             'AMY': 'NT',
             'ANA': 'KTRI',
             'ANB': 'X',
             'ANC': 'NTSI',
             'AND': 'X',
             'ANG': 'KTRM',
             'ANH': 'X',
             'ANK': 'X',
             'ANM': 'X',
             'ANN': 'X',
             'ANR': 'X',
             'ANS': 'X',
             'ANT': 'NTSI',
             'ANV': 'X',
             'ANW': 'X',
             'ANY': 'NTSI',
             'ARA': 'KR',
             'ARB': 'NKSR',
             'ARC': 'NS',
             'ARD': 'KNRS',
             'ARG': 'KR',
             'ARH': 'KNRS',
             'ARK': 'KNRS',
             'ARM': 'KNRS',
             'ARN': 'KNRS',
             'ARR': 'KR',
             'ARS': 'NKSR',
             'ART': 'NS',
             'ARV': 'KNRS',
             'ARW': 'KNRS',
             'ARY': 'NS',
             'ASA': 'TR',
             'ASB': 'TSR',
             'ASC': 'TS',
             'ASD': 'TRS',
             'ASG': 'TR',
             'ASH': 'TRS',
             'ASK': 'TRS',
             'ASM': 'TRS',
             'ASN': 'TRS',
             'ASR': 'TR',
             'ASS': 'TSR',
             'AST': 'TS',
             'ASV': 'TRS',
             'ASW': 'TRS',
             'ASY': 'TS',
             'ATA': 'I',
             'ATB': 'IM',
             'ATC': 'I',
             'ATD': 'IM',
             'ATG': 'M',
             'ATH': 'I',
             'ATK': 'MI',
             'ATM': 'I',
             'ATN': 'IM',
             'ATR': 'IM',
             'ATS': 'IM',
             'ATT': 'I',
             'ATV': 'IM',
             'ATW': 'I',
             'ATY': 'I',
             'AVA': 'KTR',
             'AVB': 'X',
             'AVC': 'NTS',
             'AVD': 'X',
             'AVG': 'KTR',
             'AVH': 'X',
             'AVK': 'X',
             'AVM': 'X',
             'AVN': 'X',
             'AVR': 'KTR',
             'AVS': 'X',
             'AVT': 'NTS',
             'AVV': 'X',
             'AVW': 'X',
             'AVY': 'NTS',
             'AWA': 'KI',
             'AWB': 'NKIM',
             'AWC': 'NI',
             'AWD': 'KNIM',
             'AWG': 'KM',
             'AWH': 'KNI',
             'AWK': 'KNMI',
             'AWM': 'KNI',
             'AWN': 'KNIM',
             'AWR': 'KIM',
             'AWS': 'NKIM',
             'AWT': 'NI',
             'AWV': 'KNIM',
             'AWW': 'KNI',
             'AWY': 'NI',
             'AYA': 'TI',
             'AYB': 'TIM',
             'AYC': 'TI',
             'AYD': 'TIM',
             'AYG': 'TM',
             'AYH': 'TI',
             'AYK': 'TMI',
             'AYM': 'TI',
             'AYN': 'TIM',
             'AYR': 'TIM',
             'AYS': 'TIM',
             'AYT': 'TI',
             'AYV': 'TIM',
             'AYW': 'TI',
             'AYY': 'TI',
             'BAA': 'QE*',
             'BAB': 'X',
             'BAC': 'HDY',
             'BAD': 'X',
             'BAG': 'QE*',
             'BAH': 'X',
             'BAK': 'X',
             'BAM': 'X',
             'BAN': 'X',
             'BAR': 'QE*',
             'BAS': 'X',
             'BAT': 'HDY',
             'BAV': 'X',
             'BAW': 'X',
             'BAY': 'HDY',
             'BBA': 'X',
             'BBB': 'X',
             'BBC': 'X',
             'BBD': 'X',
             'BBG': 'X',
             'BBH': 'X',
             'BBK': 'X',
             'BBM': 'X',
             'BBN': 'X',
             'BBR': 'X',
             'BBS': 'X',
             'BBT': 'X',
             'BBV': 'X',
             'BBW': 'X',
             'BBY': 'X',
             'BCA': 'PAS',
             'BCB': 'PAS',
             'BCC': 'PAS',
             'BCD': 'PAS',
             'BCG': 'PAS',
             'BCH': 'PAS',
             'BCK': 'PAS',
             'BCM': 'PAS',
             'BCN': 'PAS',
             'BCR': 'PAS',
             'BCS': 'PAS',
             'BCT': 'PAS',
             'BCV': 'PAS',
             'BCW': 'PAS',
             'BCY': 'PAS',
             'BDA': 'X',
             'BDB': 'X',
             'BDC': 'X',
             'BDD': 'X',
             'BDG': 'X',
             'BDH': 'X',
             'BDK': 'X',
             'BDM': 'X',
             'BDN': 'X',
             'BDR': 'X',
             'BDS': 'X',
             'BDT': 'X',
             'BDV': 'X',
             'BDW': 'X',
             'BDY': 'X',
             'BGA': 'RG*',
             'BGB': 'RGCW',
             'BGC': 'RGC',
             'BGD': 'X',
             'BGG': 'RGW',
             'BGH': 'RG*C',
             'BGK': 'RGWC',
             'BGM': 'RG*C',
             'BGN': 'X',
             'BGR': 'RG*W',
             'BGS': 'RGCW',
             'BGT': 'RGC',
             'BGV': 'X',
             'BGW': 'RG*C',
             'BGY': 'RGC',
             'BHA': 'X',
             'BHB': 'X',
             'BHC': 'X',
             'BHD': 'X',
             'BHG': 'X',
             'BHH': 'X',
             'BHK': 'X',
             'BHM': 'X',
             'BHN': 'X',
             'BHR': 'X',
             'BHS': 'X',
             'BHT': 'X',
             'BHV': 'X',
             'BHW': 'X',
             'BHY': 'X',
             'BKA': 'X',
             'BKB': 'X',
             'BKC': 'X',
             'BKD': 'X',
             'BKG': 'X',
             'BKH': 'X',
             'BKK': 'X',
             'BKM': 'X',
             'BKN': 'X',
             'BKR': 'X',
             'BKS': 'X',
             'BKT': 'X',
             'BKV': 'X',
             'BKW': 'X',
             'BKY': 'X',
             'BMA': 'X',
             'BMB': 'X',
             'BMC': 'X',
             'BMD': 'X',
             'BMG': 'X',
             'BMH': 'X',
             'BMK': 'X',
             'BMM': 'X',
             'BMN': 'X',
             'BMR': 'X',
             'BMS': 'X',
             'BMT': 'X',
             'BMV': 'X',
             'BMW': 'X',
             'BMY': 'X',
             'BNA': 'X',
             'BNB': 'X',
             'BNC': 'X',
             'BND': 'X',
             'BNG': 'X',
             'BNH': 'X',
             'BNK': 'X',
             'BNM': 'X',
             'BNN': 'X',
             'BNR': 'X',
             'BNS': 'X',
             'BNT': 'X',
             'BNV': 'X',
             'BNW': 'X',
             'BNY': 'X',
             'BRA': 'X',
             'BRB': 'X',
             'BRC': 'X',
             'BRD': 'X',
             'BRG': 'X',
             'BRH': 'X',
             'BRK': 'X',
             'BRM': 'X',
             'BRN': 'X',
             'BRR': 'X',
             'BRS': 'X',
             'BRT': 'X',
             'BRV': 'X',
             'BRW': 'X',
             'BRY': 'X',
             'BSA': 'X',
             'BSB': 'X',
             'BSC': 'X',
             'BSD': 'X',
             'BSG': 'X',
             'BSH': 'X',
             'BSK': 'X',
             'BSM': 'X',
             'BSN': 'X',
             'BSR': 'X',
             'BSS': 'X',
             'BST': 'X',
             'BSV': 'X',
             'BSW': 'X',
             'BSY': 'X',
             'BTA': 'LV',
             'BTB': 'LVF',
             'BTC': 'LVF',
             'BTD': 'LVF',
             'BTG': 'LV',
             'BTH': 'LVF',
             'BTK': 'LVF',
             'BTM': 'LVF',
             'BTN': 'LVF',
             'BTR': 'LV',
             'BTS': 'LVF',
             'BTT': 'LVF',
             'BTV': 'LVF',
             'BTW': 'LVF',
             'BTY': 'LVF',
             'BVA': 'X',
             'BVB': 'X',
             'BVC': 'X',
             'BVD': 'X',
             'BVG': 'X',
             'BVH': 'X',
             'BVK': 'X',
             'BVM': 'X',
             'BVN': 'X',
             'BVR': 'X',
             'BVS': 'X',
             'BVT': 'X',
             'BVV': 'X',
             'BVW': 'X',
             'BVY': 'X',
             'BWA': 'X',
             'BWB': 'X',
             'BWC': 'X',
             'BWD': 'X',
             'BWG': 'X',
             'BWH': 'X',
             'BWK': 'X',
             'BWM': 'X',
             'BWN': 'X',
             'BWR': 'X',
             'BWS': 'X',
             'BWT': 'X',
             'BWV': 'X',
             'BWW': 'X',
             'BWY': 'X',
             'BYA': 'X',
             'BYB': 'X',
             'BYC': 'X',
             'BYD': 'X',
             'BYG': 'X',
             'BYH': 'X',
             'BYK': 'X',
             'BYM': 'X',
             'BYN': 'X',
             'BYR': 'X',
             'BYS': 'X',
             'BYT': 'X',
             'BYV': 'X',
             'BYW': 'X',
             'BYY': 'X',
             'CAA': 'Q',
             'CAB': 'HQ',
             'CAC': 'H',
             'CAD': 'QH',
             'CAG': 'Q',
             'CAH': 'QH',
             'CAK': 'QH',
             'CAM': 'QH',
             'CAN': 'QH',
             'CAR': 'Q',
             'CAS': 'HQ',
             'CAT': 'H',
             'CAV': 'QH',
             'CAW': 'QH',
             'CAY': 'H',
             'CBA': 'PRL',
             'CBB': 'PRL',
             'CBC': 'PRL',
             'CBD': 'PRL',
             'CBG': 'PRL',
             'CBH': 'PRL',
             'CBK': 'PRL',
             'CBM': 'PRL',
             'CBN': 'PRL',
             'CBR': 'PRL',
             'CBS': 'PRL',
             'CBT': 'PRL',
             'CBV': 'PRL',
             'CBW': 'PRL',
             'CBY': 'PRL',
             'CCA': 'P',
             'CCB': 'P',
             'CCC': 'P',
             'CCD': 'P',
             'CCG': 'P',
             'CCH': 'P',
             'CCK': 'P',
             'CCM': 'P',
             'CCN': 'P',
             'CCR': 'P',
             'CCS': 'P',
             'CCT': 'P',
             'CCV': 'P',
             'CCW': 'P',
             'CCY': 'P',
             'CDA': 'QRL',
             'CDB': 'HQRL',
             'CDC': 'HRL',
             'CDD': 'QHRL',
             'CDG': 'QRL',
             'CDH': 'QHRL',
             'CDK': 'QHRL',
             'CDM': 'QHRL',
             'CDN': 'QHRL',
             'CDR': 'QRL',
             'CDS': 'HQRL',
             'CDT': 'HRL',
             'CDV': 'QHRL',
             'CDW': 'QHRL',
             'CDY': 'HRL',
             'CGA': 'R',
             'CGB': 'R',
             'CGC': 'R',
             'CGD': 'R',
             'CGG': 'R',
             'CGH': 'R',
             'CGK': 'R',
             'CGM': 'R',
             'CGN': 'R',
             'CGR': 'R',
             'CGS': 'R',
             'CGT': 'R',
             'CGV': 'R',
             'CGW': 'R',
             'CGY': 'R',
             'CHA': 'QPL',
             'CHB': 'HQPL',
             'CHC': 'HPL',
             'CHD': 'QHPL',
             'CHG': 'QPL',
             'CHH': 'QHPL',
             'CHK': 'QHPL',
             'CHM': 'QHPL',
             'CHN': 'QHPL',
             'CHR': 'QPL',
             'CHS': 'HQPL',
             'CHT': 'HPL',
             'CHV': 'QHPL',
             'CHW': 'QHPL',
             'CHY': 'HPL',
             'CKA': 'RL',
             'CKB': 'RL',
             'CKC': 'RL',
             'CKD': 'RL',
             'CKG': 'RL',
             'CKH': 'RL',
             'CKK': 'RL',
             'CKM': 'RL',
             'CKN': 'RL',
             'CKR': 'RL',
             'CKS': 'RL',
             'CKT': 'RL',
             'CKV': 'RL',
             'CKW': 'RL',
             'CKY': 'RL',
             'CMA': 'QP',
             'CMB': 'HQP',
             'CMC': 'HP',
             'CMD': 'QHP',
             'CMG': 'QP',
             'CMH': 'QHP',
             'CMK': 'QHP',
             'CMM': 'QHP',
             'CMN': 'QHP',
             'CMR': 'QP',
             'CMS': 'HQP',
             'CMT': 'HP',
             'CMV': 'QHP',
             'CMW': 'QHP',
             'CMY': 'HP',
             'CNA': 'QPRL',
             'CNB': 'X',
             'CNC': 'HPRL',
             'CND': 'X',
             'CNG': 'QPRL',
             'CNH': 'X',
             'CNK': 'X',
             'CNM': 'X',
             'CNN': 'X',
             'CNR': 'QPRL',
             'CNS': 'X',
             'CNT': 'HPRL',
             'CNV': 'X',
             'CNW': 'X',
             'CNY': 'HPRL',
             'CRA': 'QR',
             'CRB': 'HQR',
             'CRC': 'HR',
             'CRD': 'QHR',
             'CRG': 'QR',
             'CRH': 'QHR',
             'CRK': 'QHR',
             'CRM': 'QHR',
             'CRN': 'QHR',
             'CRR': 'QR',
             'CRS': 'HQR',
             'CRT': 'HR',
             'CRV': 'QHR',
             'CRW': 'QHR',
             'CRY': 'HR',
             'CSA': 'PR',
             'CSB': 'PR',
             'CSC': 'PR',
             'CSD': 'PR',
             'CSG': 'PR',
             'CSH': 'PR',
             'CSK': 'PR',
             'CSM': 'PR',
             'CSN': 'PR',
             'CSR': 'PR',
             'CSS': 'PR',
             'CST': 'PR',
             'CSV': 'PR',
             'CSW': 'PR',
             'CSY': 'PR',
             'CTA': 'L',
             'CTB': 'L',
             'CTC': 'L',
             'CTD': 'L',
             'CTG': 'L',
             'CTH': 'L',
             'CTK': 'L',
             'CTM': 'L',
             'CTN': 'L',
             'CTR': 'L',
             'CTS': 'L',
             'CTT': 'L',
             'CTV': 'L',
             'CTW': 'L',
             'CTY': 'L',
             'CVA': 'QPR',
             'CVB': 'HQPR',
             'CVC': 'HPR',
             'CVD': 'QHPR',
             'CVG': 'QPR',
             'CVH': 'QHPR',
             'CVK': 'QHPR',
             'CVM': 'QHPR',
             'CVN': 'QHPR',
             'CVR': 'QPR',
             'CVS': 'HQPR',
             'CVT': 'HPR',
             'CVV': 'QHPR',
             'CVW': 'QHPR',
             'CVY': 'HPR',
             'CWA': 'QL',
             'CWB': 'HQL',
             'CWC': 'HL',
             'CWD': 'QHL',
             'CWG': 'QL',
             'CWH': 'QHL',
             'CWK': 'QHL',
             'CWM': 'QHL',
             'CWN': 'QHL',
             'CWR': 'QL',
             'CWS': 'HQL',
             'CWT': 'HL',
             'CWV': 'QHL',
             'CWW': 'QHL',
             'CWY': 'HL',
             'CYA': 'PL',
             'CYB': 'PL',
             'CYC': 'PL',
             'CYD': 'PL',
             'CYG': 'PL',
             'CYH': 'PL',
             'CYK': 'PL',
             'CYM': 'PL',
             'CYN': 'PL',
             'CYR': 'PL',
             'CYS': 'PL',
             'CYT': 'PL',
             'CYV': 'PL',
             'CYW': 'PL',
             'CYY': 'PL',
             'DAA': 'KE*',
             'DAB': 'X',
             'DAC': 'NDY',
             'DAD': 'X',
             'DAG': 'KE*',
             'DAH': 'X',
             'DAK': 'X',
             'DAM': 'X',
             'DAN': 'X',
             'DAR': 'KE*',
             'DAS': 'X',
             'DAT': 'NDY',
             'DAV': 'X',
             'DAW': 'X',
             'DAY': 'NDY',
             'DBA': 'X',
             'DBB': 'X',
             'DBC': 'X',
             'DBD': 'X',
             'DBG': 'X',
             'DBH': 'X',
             'DBK': 'X',
             'DBM': 'X',
             'DBN': 'X',
             'DBR': 'X',
             'DBS': 'X',
             'DBT': 'X',
             'DBV': 'X',
             'DBW': 'X',
             'DBY': 'X',
             'DCA': 'TAS',
             'DCB': 'TAS',
             'DCC': 'TAS',
             'DCD': 'TAS',
             'DCG': 'TAS',
             'DCH': 'TAS',
             'DCK': 'TAS',
             'DCM': 'TAS',
             'DCN': 'TAS',
             'DCR': 'TAS',
             'DCS': 'TAS',
             'DCT': 'TAS',
             'DCV': 'TAS',
             'DCW': 'TAS',
             'DCY': 'TAS',
             'DDA': 'X',
             'DDB': 'X',
             'DDC': 'X',
             'DDD': 'X',
             'DDG': 'X',
             'DDH': 'X',
             'DDK': 'X',
             'DDM': 'X',
             'DDN': 'X',
             'DDR': 'X',
             'DDS': 'X',
             'DDT': 'X',
             'DDV': 'X',
             'DDW': 'X',
             'DDY': 'X',
             'DGA': 'RG*',
             'DGB': 'X',
             'DGC': 'SGC',
             'DGD': 'X',
             'DGG': 'RGW',
             'DGH': 'X',
             'DGK': 'X',
             'DGM': 'X',
             'DGN': 'X',
             'DGR': 'RG*W',
             'DGS': 'X',
             'DGT': 'SGC',
             'DGV': 'X',
             'DGW': 'X',
             'DGY': 'SGC',
             'DHA': 'X',
             'DHB': 'X',
             'DHC': 'X',
             'DHD': 'X',
             'DHG': 'X',
             'DHH': 'X',
             'DHK': 'X',
             'DHM': 'X',
             'DHN': 'X',
             'DHR': 'X',
             'DHS': 'X',
             'DHT': 'X',
             'DHV': 'X',
             'DHW': 'X',
             'DHY': 'X',
             'DKA': 'X',
             'DKB': 'X',
             'DKC': 'X',
             'DKD': 'X',
             'DKG': 'X',
             'DKH': 'X',
             'DKK': 'X',
             'DKM': 'X',
             'DKN': 'X',
             'DKR': 'X',
             'DKS': 'X',
             'DKT': 'X',
             'DKV': 'X',
             'DKW': 'X',
             'DKY': 'X',
             'DMA': 'X',
             'DMB': 'X',
             'DMC': 'X',
             'DMD': 'X',
             'DMG': 'X',
             'DMH': 'X',
             'DMK': 'X',
             'DMM': 'X',
             'DMN': 'X',
             'DMR': 'X',
             'DMS': 'X',
             'DMT': 'X',
             'DMV': 'X',
             'DMW': 'X',
             'DMY': 'X',
             'DNA': 'X',
             'DNB': 'X',
             'DNC': 'X',
             'DND': 'X',
             'DNG': 'X',
             'DNH': 'X',
             'DNK': 'X',
             'DNM': 'X',
             'DNN': 'X',
             'DNR': 'X',
             'DNS': 'X',
             'DNT': 'X',
             'DNV': 'X',
             'DNW': 'X',
             'DNY': 'X',
             'DRA': 'X',
             'DRB': 'X',
             'DRC': 'X',
             'DRD': 'X',
             'DRG': 'X',
             'DRH': 'X',
             'DRK': 'X',
             'DRM': 'X',
             'DRN': 'X',
             'DRR': 'X',
             'DRS': 'X',
             'DRT': 'X',
             'DRV': 'X',
             'DRW': 'X',
             'DRY': 'X',
             'DSA': 'X',
             'DSB': 'X',
             'DSC': 'X',
             'DSD': 'X',
             'DSG': 'X',
             'DSH': 'X',
             'DSK': 'X',
             'DSM': 'X',
             'DSN': 'X',
             'DSR': 'X',
             'DSS': 'X',
             'DST': 'X',
             'DSV': 'X',
             'DSW': 'X',
             'DSY': 'X',
             'DTA': 'IVL',
             'DTB': 'X',
             'DTC': 'IVF',
             'DTD': 'X',
             'DTG': 'MVL',
             'DTH': 'IVLF',
             'DTK': 'X',
             'DTM': 'IVLF',
             'DTN': 'X',
             'DTR': 'IMVL',
             'DTS': 'X',
             'DTT': 'IVF',
             'DTV': 'X',
             'DTW': 'IVLF',
             'DTY': 'IVF',
             'DVA': 'X',
             'DVB': 'X',
             'DVC': 'X',
             'DVD': 'X',
             'DVG': 'X',
             'DVH': 'X',
             'DVK': 'X',
             'DVM': 'X',
             'DVN': 'X',
             'DVR': 'X',
             'DVS': 'X',
             'DVT': 'X',
             'DVV': 'X',
             'DVW': 'X',
             'DVY': 'X',
             'DWA': 'X',
             'DWB': 'X',
             'DWC': 'X',
             'DWD': 'X',
             'DWG': 'X',
             'DWH': 'X',
             'DWK': 'X',
             'DWM': 'X',
             'DWN': 'X',
             'DWR': 'X',
             'DWS': 'X',
             'DWT': 'X',
             'DWV': 'X',
             'DWW': 'X',
             'DWY': 'X',
             'DYA': 'X',
             'DYB': 'X',
             'DYC': 'X',
             'DYD': 'X',
             'DYG': 'X',
             'DYH': 'X',
             'DYK': 'X',
             'DYM': 'X',
             'DYN': 'X',
             'DYR': 'X',
             'DYS': 'X',
             'DYT': 'X',
             'DYV': 'X',
             'DYW': 'X',
             'DYY': 'X',
             'GAA': 'E',
             'GAB': 'DE',
             'GAC': 'D',
             'GAD': 'ED',
             'GAG': 'E',
             'GAH': 'ED',
             'GAK': 'ED',
             'GAM': 'ED',
             'GAN': 'ED',
             'GAR': 'E',
             'GAS': 'DE',
             'GAT': 'D',
             'GAV': 'ED',
             'GAW': 'ED',
             'GAY': 'D',
             'GBA': 'AGV',
             'GBB': 'AGV',
             'GBC': 'AGV',
             'GBD': 'AGV',
             'GBG': 'AGV',
             'GBH': 'AGV',
             'GBK': 'AGV',
             'GBM': 'AGV',
             'GBN': 'AGV',
             'GBR': 'AGV',
             'GBS': 'AGV',
             'GBT': 'AGV',
             'GBV': 'AGV',
             'GBW': 'AGV',
             'GBY': 'AGV',
             'GCA': 'A',
             'GCB': 'A',
             'GCC': 'A',
             'GCD': 'A',
             'GCG': 'A',
             'GCH': 'A',
             'GCK': 'A',
             'GCM': 'A',
             'GCN': 'A',
             'GCR': 'A',
             'GCS': 'A',
             'GCT': 'A',
             'GCV': 'A',
             'GCW': 'A',
             'GCY': 'A',
             'GDA': 'EGV',
             'GDB': 'DEGV',
             'GDC': 'DGV',
             'GDD': 'EDGV',
             'GDG': 'EGV',
             'GDH': 'EDGV',
             'GDK': 'EDGV',
             'GDM': 'EDGV',
             'GDN': 'EDGV',
             'GDR': 'EGV',
             'GDS': 'DEGV',
             'GDT': 'DGV',
             'GDV': 'EDGV',
             'GDW': 'EDGV',
             'GDY': 'DGV',
             'GGA': 'G',
             'GGB': 'G',
             'GGC': 'G',
             'GGD': 'G',
             'GGG': 'G',
             'GGH': 'G',
             'GGK': 'G',
             'GGM': 'G',
             'GGN': 'G',
             'GGR': 'G',
             'GGS': 'G',
             'GGT': 'G',
             'GGV': 'G',
             'GGW': 'G',
             'GGY': 'G',
             'GHA': 'EAV',
             'GHB': 'DEAV',
             'GHC': 'DAV',
             'GHD': 'EDAV',
             'GHG': 'EAV',
             'GHH': 'EDAV',
             'GHK': 'EDAV',
             'GHM': 'EDAV',
             'GHN': 'EDAV',
             'GHR': 'EAV',
             'GHS': 'DEAV',
             'GHT': 'DAV',
             'GHV': 'EDAV',
             'GHW': 'EDAV',
             'GHY': 'DAV',
             'GKA': 'GV',
             'GKB': 'GV',
             'GKC': 'GV',
             'GKD': 'GV',
             'GKG': 'GV',
             'GKH': 'GV',
             'GKK': 'GV',
             'GKM': 'GV',
             'GKN': 'GV',
             'GKR': 'GV',
             'GKS': 'GV',
             'GKT': 'GV',
             'GKV': 'GV',
             'GKW': 'GV',
             'GKY': 'GV',
             'GMA': 'EA',
             'GMB': 'DEA',
             'GMC': 'DA',
             'GMD': 'EDA',
             'GMG': 'EA',
             'GMH': 'EDA',
             'GMK': 'EDA',
             'GMM': 'EDA',
             'GMN': 'EDA',
             'GMR': 'EA',
             'GMS': 'DEA',
             'GMT': 'DA',
             'GMV': 'EDA',
             'GMW': 'EDA',
             'GMY': 'DA',
             'GNA': 'EAGV',
             'GNB': 'X',
             'GNC': 'DAGV',
             'GND': 'X',
             'GNG': 'EAGV',
             'GNH': 'X',
             'GNK': 'X',
             'GNM': 'X',
             'GNN': 'X',
             'GNR': 'EAGV',
             'GNS': 'X',
             'GNT': 'DAGV',
             'GNV': 'X',
             'GNW': 'X',
             'GNY': 'DAGV',
             'GRA': 'EG',
             'GRB': 'DEG',
             'GRC': 'DG',
             'GRD': 'EDG',
             'GRG': 'EG',
             'GRH': 'EDG',
             'GRK': 'EDG',
             'GRM': 'EDG',
             'GRN': 'EDG',
             'GRR': 'EG',
             'GRS': 'DEG',
             'GRT': 'DG',
             'GRV': 'EDG',
             'GRW': 'EDG',
             'GRY': 'DG',
             'GSA': 'AG',
             'GSB': 'AG',
             'GSC': 'AG',
             'GSD': 'AG',
             'GSG': 'AG',
             'GSH': 'AG',
             'GSK': 'AG',
             'GSM': 'AG',
             'GSN': 'AG',
             'GSR': 'AG',
             'GSS': 'AG',
             'GST': 'AG',
             'GSV': 'AG',
             'GSW': 'AG',
             'GSY': 'AG',
             'GTA': 'V',
             'GTB': 'V',
             'GTC': 'V',
             'GTD': 'V',
             'GTG': 'V',
             'GTH': 'V',
             'GTK': 'V',
             'GTM': 'V',
             'GTN': 'V',
             'GTR': 'V',
             'GTS': 'V',
             'GTT': 'V',
             'GTV': 'V',
             'GTW': 'V',
             'GTY': 'V',
             'GVA': 'EAG',
             'GVB': 'DEAG',
             'GVC': 'DAG',
             'GVD': 'EDAG',
             'GVG': 'EAG',
             'GVH': 'EDAG',
             'GVK': 'EDAG',
             'GVM': 'EDAG',
             'GVN': 'EDAG',
             'GVR': 'EAG',
             'GVS': 'DEAG',
             'GVT': 'DAG',
             'GVV': 'EDAG',
             'GVW': 'EDAG',
             'GVY': 'DAG',
             'GWA': 'EV',
             'GWB': 'DEV',
             'GWC': 'DV',
             'GWD': 'EDV',
             'GWG': 'EV',
             'GWH': 'EDV',
             'GWK': 'EDV',
             'GWM': 'EDV',
             'GWN': 'EDV',
             'GWR': 'EV',
             'GWS': 'DEV',
             'GWT': 'DV',
             'GWV': 'EDV',
             'GWW': 'EDV',
             'GWY': 'DV',
             'GYA': 'AV',
             'GYB': 'AV',
             'GYC': 'AV',
             'GYD': 'AV',
             'GYG': 'AV',
             'GYH': 'AV',
             'GYK': 'AV',
             'GYM': 'AV',
             'GYN': 'AV',
             'GYR': 'AV',
             'GYS': 'AV',
             'GYT': 'AV',
             'GYV': 'AV',
             'GYW': 'AV',
             'GYY': 'AV',
             'HAA': 'KQ*',
             'HAB': 'X',
             'HAC': 'NHY',
             'HAD': 'X',
             'HAG': 'KQ*',
             'HAH': 'X',
             'HAK': 'X',
             'HAM': 'X',
             'HAN': 'X',
             'HAR': 'KQ*',
             'HAS': 'X',
             'HAT': 'NHY',
             'HAV': 'X',
             'HAW': 'X',
             'HAY': 'NHY',
             'HBA': 'X',
             'HBB': 'X',
             'HBC': 'X',
             'HBD': 'X',
             'HBG': 'X',
             'HBH': 'X',
             'HBK': 'X',
             'HBM': 'X',
             'HBN': 'X',
             'HBR': 'X',
             'HBS': 'X',
             'HBT': 'X',
             'HBV': 'X',
             'HBW': 'X',
             'HBY': 'X',
             'HCA': 'TPS',
             'HCB': 'TPS',
             'HCC': 'TPS',
             'HCD': 'TPS',
             'HCG': 'TPS',
             'HCH': 'TPS',
             'HCK': 'TPS',
             'HCM': 'TPS',
             'HCN': 'TPS',
             'HCR': 'TPS',
             'HCS': 'TPS',
             'HCT': 'TPS',
             'HCV': 'TPS',
             'HCW': 'TPS',
             'HCY': 'TPS',
             'HDA': 'X',
             'HDB': 'X',
             'HDC': 'X',
             'HDD': 'X',
             'HDG': 'X',
             'HDH': 'X',
             'HDK': 'X',
             'HDM': 'X',
             'HDN': 'X',
             'HDR': 'X',
             'HDS': 'X',
             'HDT': 'X',
             'HDV': 'X',
             'HDW': 'X',
             'HDY': 'X',
             'HGA': 'R*',
             'HGB': 'SRCW',
             'HGC': 'SRC',
             'HGD': 'X',
             'HGG': 'RW',
             'HGH': 'RS*C',
             'HGK': 'RSWC',
             'HGM': 'RS*C',
             'HGN': 'X',
             'HGR': 'R*W',
             'HGS': 'SRCW',
             'HGT': 'SRC',
             'HGV': 'X',
             'HGW': 'RS*C',
             'HGY': 'SRC',
             'HHA': 'X',
             'HHB': 'X',
             'HHC': 'X',
             'HHD': 'X',
             'HHG': 'X',
             'HHH': 'X',
             'HHK': 'X',
             'HHM': 'X',
             'HHN': 'X',
             'HHR': 'X',
             'HHS': 'X',
             'HHT': 'X',
             'HHV': 'X',
             'HHW': 'X',
             'HHY': 'X',
             'HKA': 'RIL*',
             'HKB': 'X',
             'HKC': 'X',
             'HKD': 'X',
             'HKG': 'RMLW',
             'HKH': 'X',
             'HKK': 'X',
             'HKM': 'X',
             'HKN': 'X',
             'HKR': 'X',
             'HKS': 'X',
             'HKT': 'X',
             'HKV': 'X',
             'HKW': 'X',
             'HKY': 'X',
             'HMA': 'X',
             'HMB': 'X',
             'HMC': 'X',
             'HMD': 'X',
             'HMG': 'X',
             'HMH': 'X',
             'HMK': 'X',
             'HMM': 'X',
             'HMN': 'X',
             'HMR': 'X',
             'HMS': 'X',
             'HMT': 'X',
             'HMV': 'X',
             'HMW': 'X',
             'HMY': 'X',
             'HNA': 'X',
             'HNB': 'X',
             'HNC': 'X',
             'HND': 'X',
             'HNG': 'X',
             'HNH': 'X',
             'HNK': 'X',
             'HNM': 'X',
             'HNN': 'X',
             'HNR': 'X',
             'HNS': 'X',
             'HNT': 'X',
             'HNV': 'X',
             'HNW': 'X',
             'HNY': 'X',
             'HRA': 'KRQ*',
             'HRB': 'X',
             'HRC': 'X',
             'HRD': 'X',
             'HRG': 'X',
             'HRH': 'X',
             'HRK': 'X',
             'HRM': 'X',
             'HRN': 'X',
             'HRR': 'X',
             'HRS': 'X',
             'HRT': 'X',
             'HRV': 'X',
             'HRW': 'X',
             'HRY': 'X',
             'HSA': 'X',
             'HSB': 'X',
             'HSC': 'X',
             'HSD': 'X',
             'HSG': 'X',
             'HSH': 'X',
             'HSK': 'X',
             'HSM': 'X',
             'HSN': 'X',
             'HSR': 'X',
             'HSS': 'X',
             'HST': 'X',
             'HSV': 'X',
             'HSW': 'X',
             'HSY': 'X',
             'HTA': 'IL',
             'HTB': 'IMLF',
             'HTC': 'ILF',
             'HTD': 'IMLF',
             'HTG': 'ML',
             'HTH': 'ILF',
             'HTK': 'MILF',
             'HTM': 'ILF',
             'HTN': 'IMLF',
             'HTR': 'IML',
             'HTS': 'IMLF',
             'HTT': 'ILF',
             'HTV': 'IMLF',
             'HTW': 'ILF',
             'HTY': 'ILF',
             'HVA': 'X',
             'HVB': 'X',
             'HVC': 'X',
             'HVD': 'X',
             'HVG': 'X',
             'HVH': 'X',
             'HVK': 'X',
             'HVM': 'X',
             'HVN': 'X',
             'HVR': 'X',
             'HVS': 'X',
             'HVT': 'X',
             'HVV': 'X',
             'HVW': 'X',
             'HVY': 'X',
             'HWA': 'X',
             'HWB': 'X',
             'HWC': 'X',
             'HWD': 'X',
             'HWG': 'X',
             'HWH': 'X',
             'HWK': 'X',
             'HWM': 'X',
             'HWN': 'X',
             'HWR': 'X',
             'HWS': 'X',
             'HWT': 'X',
             'HWV': 'X',
             'HWW': 'X',
             'HWY': 'X',
             'HYA': 'X',
             'HYB': 'X',
             'HYC': 'X',
             'HYD': 'X',
             'HYG': 'X',
             'HYH': 'X',
             'HYK': 'X',
             'HYM': 'X',
             'HYN': 'X',
             'HYR': 'X',
             'HYS': 'X',
             'HYT': 'X',
             'HYV': 'X',
             'HYW': 'X',
             'HYY': 'X',
             'KAA': 'E*',
             'KAB': 'DEY*',
             'KAC': 'DY',
             'KAD': 'ED*Y',
             'KAG': 'E*',
             'KAH': 'ED*Y',
             'KAK': 'ED*Y',
             'KAM': 'ED*Y',
             'KAN': 'ED*Y',
             'KAR': 'E*',
             'KAS': 'DEY*',
             'KAT': 'DY',
             'KAV': 'ED*Y',
             'KAW': 'ED*Y',
             'KAY': 'DY',
             'KBA': 'X',
             'KBB': 'X',
             'KBC': 'X',
             'KBD': 'X',
             'KBG': 'X',
             'KBH': 'X',
             'KBK': 'X',
             'KBM': 'X',
             'KBN': 'X',
             'KBR': 'X',
             'KBS': 'X',
             'KBT': 'X',
             'KBV': 'X',
             'KBW': 'X',
             'KBY': 'X',
             'KCA': 'AS',
             'KCB': 'AS',
             'KCC': 'AS',
             'KCD': 'AS',
             'KCG': 'AS',
             'KCH': 'AS',
             'KCK': 'AS',
             'KCM': 'AS',
             'KCN': 'AS',
             'KCR': 'AS',
             'KCS': 'AS',
             'KCT': 'AS',
             'KCV': 'AS',
             'KCW': 'AS',
             'KCY': 'AS',
             'KDA': 'X',
             'KDB': 'X',
             'KDC': 'X',
             'KDD': 'X',
             'KDG': 'X',
             'KDH': 'X',
             'KDK': 'X',
             'KDM': 'X',
             'KDN': 'X',
             'KDR': 'X',
             'KDS': 'X',
             'KDT': 'X',
             'KDV': 'X',
             'KDW': 'X',
             'KDY': 'X',
             'KGA': 'G*',
             'KGB': 'GCW',
             'KGC': 'GC',
             'KGD': 'G*WC',
             'KGG': 'GW',
             'KGH': 'G*C',
             'KGK': 'GWC',
             'KGM': 'G*C',
             'KGN': 'G*CW',
             'KGR': 'G*W',
             'KGS': 'GCW',
             'KGT': 'GC',
             'KGV': 'G*CW',
             'KGW': 'G*C',
             'KGY': 'GC',
             'KHA': 'X',
             'KHB': 'X',
             'KHC': 'X',
             'KHD': 'X',
             'KHG': 'X',
             'KHH': 'X',
             'KHK': 'X',
             'KHM': 'X',
             'KHN': 'X',
             'KHR': 'X',
             'KHS': 'X',
             'KHT': 'X',
             'KHV': 'X',
             'KHW': 'X',
             'KHY': 'X',
             'KKA': 'GV*L',
             'KKB': 'X',
             'KKC': 'GVCF',
             'KKD': 'X',
             'KKG': 'GVWL',
             'KKH': 'X',
             'KKK': 'X',
             'KKM': 'X',
             'KKN': 'X',
             'KKR': 'X',
             'KKS': 'X',
             'KKT': 'GVCF',
             'KKV': 'X',
             'KKW': 'X',
             'KKY': 'GVCF',
             'KMA': 'EA*S',
             'KMB': 'X',
             'KMC': 'DAYS',
             'KMD': 'X',
             'KMG': 'EA*S',
             'KMH': 'X',
             'KMK': 'X',
             'KMM': 'X',
             'KMN': 'X',
             'KMR': 'EA*S',
             'KMS': 'X',
             'KMT': 'DAYS',
             'KMV': 'X',
             'KMW': 'X',
             'KMY': 'DAYS',
             'KNA': 'X',
             'KNB': 'X',
             'KNC': 'X',
             'KND': 'X',
             'KNG': 'X',
             'KNH': 'X',
             'KNK': 'X',
             'KNM': 'X',
             'KNN': 'X',
             'KNR': 'X',
             'KNS': 'X',
             'KNT': 'X',
             'KNV': 'X',
             'KNW': 'X',
             'KNY': 'X',
             'KRA': 'EG*',
             'KRB': 'X',
             'KRC': 'DGYC',
             'KRD': 'X',
             'KRG': 'EG*W',
             'KRH': 'X',
             'KRK': 'X',
             'KRM': 'X',
             'KRN': 'X',
             'KRR': 'EG*W',
             'KRS': 'X',
             'KRT': 'DGYC',
             'KRV': 'X',
             'KRW': 'X',
             'KRY': 'DGYC',
             'KSA': 'AGS*',
             'KSB': 'X',
             'KSC': 'AGSC',
             'KSD': 'X',
             'KSG': 'AGSW',
             'KSH': 'X',
             'KSK': 'X',
             'KSM': 'X',
             'KSN': 'X',
             'KSR': 'X',
             'KSS': 'X',
             'KST': 'AGSC',
             'KSV': 'X',
             'KSW': 'X',
             'KSY': 'AGSC',
             'KTA': 'VL',
             'KTB': 'VFL',
             'KTC': 'VF',
             'KTD': 'VLF',
             'KTG': 'VL',
             'KTH': 'VLF',
             'KTK': 'VLF',
             'KTM': 'VLF',
             'KTN': 'VLF',
             'KTR': 'VL',
             'KTS': 'VFL',
             'KTT': 'VF',
             'KTV': 'VLF',
             'KTW': 'VLF',
             'KTY': 'VF',
             'KVA': 'X',
             'KVB': 'X',
             'KVC': 'X',
             'KVD': 'X',
             'KVG': 'X',
             'KVH': 'X',
             'KVK': 'X',
             'KVM': 'X',
             'KVN': 'X',
             'KVR': 'X',
             'KVS': 'X',
             'KVT': 'X',
             'KVV': 'X',
             'KVW': 'X',
             'KVY': 'X',
             'KWA': 'EV*L',
             'KWB': 'X',
             'KWC': 'DVYF',
             'KWD': 'X',
             'KWG': 'EV*L',
             'KWH': 'X',
             'KWK': 'X',
             'KWM': 'X',
             'KWN': 'X',
             'KWR': 'EV*L',
             'KWS': 'X',
             'KWT': 'DVYF',
             'KWV': 'X',
             'KWW': 'X',
             'KWY': 'DVYF',
             'KYA': 'AVSL',
             'KYB': 'X',
             'KYC': 'AVSF',
             'KYD': 'X',
             'KYG': 'AVSL',
             'KYH': 'X',
             'KYK': 'X',
             'KYM': 'X',
             'KYN': 'X',
             'KYR': 'AVSL',
             'KYS': 'X',
             'KYT': 'AVSF',
             'KYV': 'X',
             'KYW': 'X',
             'KYY': 'AVSF',
             'MAA': 'KQ',
             'MAB': 'NKHQ',
             'MAC': 'NH',
             'MAD': 'KNQH',
             'MAG': 'KQ',
             'MAH': 'KNQH',
             'MAK': 'KNQH',
             'MAM': 'KNQH',
             'MAN': 'KNQH',
             'MAR': 'KQ',
             'MAS': 'NKHQ',
             'MAT': 'NH',
             'MAV': 'KNQH',
             'MAW': 'KNQH',
             'MAY': 'NH',
             'MBA': 'X',
             'MBB': 'X',
             'MBC': 'X',
             'MBD': 'X',
             'MBG': 'X',
             'MBH': 'X',
             'MBK': 'X',
             'MBM': 'X',
             'MBN': 'X',
             'MBR': 'X',
             'MBS': 'X',
             'MBT': 'X',
             'MBV': 'X',
             'MBW': 'X',
             'MBY': 'X',
             'MCA': 'TP',
             'MCB': 'TP',
             'MCC': 'TP',
             'MCD': 'TP',
             'MCG': 'TP',
             'MCH': 'TP',
             'MCK': 'TP',
             'MCM': 'TP',
             'MCN': 'TP',
             'MCR': 'TP',
             'MCS': 'TP',
             'MCT': 'TP',
             'MCV': 'TP',
             'MCW': 'TP',
             'MCY': 'TP',
             'MDA': 'X',
             'MDB': 'X',
             'MDC': 'X',
             'MDD': 'X',
             'MDG': 'X',
             'MDH': 'X',
             'MDK': 'X',
             'MDM': 'X',
             'MDN': 'X',
             'MDR': 'X',
             'MDS': 'X',
             'MDT': 'X',
             'MDV': 'X',
             'MDW': 'X',
             'MDY': 'X',
             'MGA': 'R',
             'MGB': 'SR',
             'MGC': 'SR',
             'MGD': 'RS',
             'MGG': 'R',
             'MGH': 'RS',
             'MGK': 'RS',
             'MGM': 'RS',
             'MGN': 'RS',
             'MGR': 'R',
             'MGS': 'SR',
             'MGT': 'SR',
             'MGV': 'RS',
             'MGW': 'RS',
             'MGY': 'SR',
             'MHA': 'X',
             'MHB': 'X',
             'MHC': 'X',
             'MHD': 'X',
             'MHG': 'X',
             'MHH': 'X',
             'MHK': 'X',
             'MHM': 'X',
             'MHN': 'X',
             'MHR': 'X',
             'MHS': 'X',
             'MHT': 'X',
             'MHV': 'X',
             'MHW': 'X',
             'MHY': 'X',
             'MKA': 'RIL',
             'MKB': 'X',
             'MKC': 'SIRL',
             'MKD': 'X',
             'MKG': 'RML',
             'MKH': 'RSIL',
             'MKK': 'X',
             'MKM': 'RSIL',
             'MKN': 'X',
             'MKR': 'RIML',
             'MKS': 'X',
             'MKT': 'SIRL',
             'MKV': 'X',
             'MKW': 'RSIL',
             'MKY': 'SIRL',
             'MMA': 'KTQP',
             'MMB': 'X',
             'MMC': 'NTHP',
             'MMD': 'X',
             'MMG': 'KTQP',
             'MMH': 'X',
             'MMK': 'X',
             'MMM': 'X',
             'MMN': 'X',
             'MMR': 'KTQP',
             'MMS': 'X',
             'MMT': 'NTHP',
             'MMV': 'X',
             'MMW': 'X',
             'MMY': 'NTHP',
             'MNA': 'X',
             'MNB': 'X',
             'MNC': 'X',
             'MND': 'X',
             'MNG': 'X',
             'MNH': 'X',
             'MNK': 'X',
             'MNM': 'X',
             'MNN': 'X',
             'MNR': 'X',
             'MNS': 'X',
             'MNT': 'X',
             'MNV': 'X',
             'MNW': 'X',
             'MNY': 'X',
             'MRA': 'KRQ',
             'MRB': 'X',
             'MRC': 'NSHR',
             'MRD': 'X',
             'MRG': 'KRQ',
             'MRH': 'X',
             'MRK': 'X',
             'MRM': 'X',
             'MRN': 'X',
             'MRR': 'KRQ',
             'MRS': 'X',
             'MRT': 'NSHR',
             'MRV': 'X',
             'MRW': 'X',
             'MRY': 'NSHR',
             'MSA': 'TRP',
             'MSB': 'TSRP',
             'MSC': 'TSPR',
             'MSD': 'TRSP',
             'MSG': 'TRP',
             'MSH': 'TRSP',
             'MSK': 'TRSP',
             'MSM': 'TRSP',
             'MSN': 'TRSP',
             'MSR': 'TRP',
             'MSS': 'TSRP',
             'MST': 'TSPR',
             'MSV': 'TRSP',
             'MSW': 'TRSP',
             'MSY': 'TSPR',
             'MTA': 'IL',
             'MTB': 'IML',
             'MTC': 'IL',
             'MTD': 'IML',
             'MTG': 'ML',
             'MTH': 'IL',
             'MTK': 'MIL',
             'MTM': 'IL',
             'MTN': 'IML',
             'MTR': 'IML',
             'MTS': 'IML',
             'MTT': 'IL',
             'MTV': 'IML',
             'MTW': 'IL',
             'MTY': 'IL',
             'MVA': 'X',
             'MVB': 'X',
             'MVC': 'X',
             'MVD': 'X',
             'MVG': 'X',
             'MVH': 'X',
             'MVK': 'X',
             'MVM': 'X',
             'MVN': 'X',
             'MVR': 'X',
             'MVS': 'X',
             'MVT': 'X',
             'MVV': 'X',
             'MVW': 'X',
             'MVY': 'X',
             'MWA': 'KIQL',
             'MWB': 'X',
             'MWC': 'NIHL',
             'MWD': 'X',
             'MWG': 'KMQL',
             'MWH': 'X',
             'MWK': 'X',
             'MWM': 'X',
             'MWN': 'X',
             'MWR': 'X',
             'MWS': 'X',
             'MWT': 'NIHL',
             'MWV': 'X',
             'MWW': 'X',
             'MWY': 'NIHL',
             'MYA': 'TIPL',
             'MYB': 'X',
             'MYC': 'TIPL',
             'MYD': 'X',
             'MYG': 'TMPL',
             'MYH': 'TIPL',
             'MYK': 'X',
             'MYM': 'TIPL',
             'MYN': 'X',
             'MYR': 'X',
             'MYS': 'X',
             'MYT': 'TIPL',
             'MYV': 'X',
             'MYW': 'TIPL',
             'MYY': 'TIPL',
             'NAA': 'KQE*',
             'NAB': 'X',
             'NAC': 'NHDY',
             'NAD': 'X',
             'NAG': 'KQE*',
             'NAH': 'X',
             'NAK': 'X',
             'NAM': 'X',
             'NAN': 'X',
             'NAR': 'KQE*',
             'NAS': 'X',
             'NAT': 'NHDY',
             'NAV': 'X',
             'NAW': 'X',
             'NAY': 'NHDY',
             'NBA': 'X',
             'NBB': 'X',
             'NBC': 'X',
             'NBD': 'X',
             'NBG': 'X',
             'NBH': 'X',
             'NBK': 'X',
             'NBM': 'X',
             'NBN': 'X',
             'NBR': 'X',
             'NBS': 'X',
             'NBT': 'X',
             'NBV': 'X',
             'NBW': 'X',
             'NBY': 'X',
             'NCA': 'TPAS',
             'NCB': 'TPAS',
             'NCC': 'TPAS',
             'NCD': 'TPAS',
             'NCG': 'TPAS',
             'NCH': 'TPAS',
             'NCK': 'TPAS',
             'NCM': 'TPAS',
             'NCN': 'TPAS',
             'NCR': 'TPAS',
             'NCS': 'TPAS',
             'NCT': 'TPAS',
             'NCV': 'TPAS',
             'NCW': 'TPAS',
             'NCY': 'TPAS',
             'NDA': 'X',
             'NDB': 'X',
             'NDC': 'X',
             'NDD': 'X',
             'NDG': 'X',
             'NDH': 'X',
             'NDK': 'X',
             'NDM': 'X',
             'NDN': 'X',
             'NDR': 'X',
             'NDS': 'X',
             'NDT': 'X',
             'NDV': 'X',
             'NDW': 'X',
             'NDY': 'X',
             'NGA': 'RG*',
             'NGB': 'X',
             'NGC': 'SRGC',
             'NGD': 'X',
             'NGG': 'RGW',
             'NGH': 'X',
             'NGK': 'X',
             'NGM': 'X',
             'NGN': 'X',
             'NGR': 'RG*W',
             'NGS': 'X',
             'NGT': 'SRGC',
             'NGV': 'X',
             'NGW': 'X',
             'NGY': 'SRGC',
             'NHA': 'X',
             'NHB': 'X',
             'NHC': 'X',
             'NHD': 'X',
             'NHG': 'X',
             'NHH': 'X',
             'NHK': 'X',
             'NHM': 'X',
             'NHN': 'X',
             'NHR': 'X',
             'NHS': 'X',
             'NHT': 'X',
             'NHV': 'X',
             'NHW': 'X',
             'NHY': 'X',
             'NKA': 'X',
             'NKB': 'X',
             'NKC': 'X',
             'NKD': 'X',
             'NKG': 'X',
             'NKH': 'X',
             'NKK': 'X',
             'NKM': 'X',
             'NKN': 'X',
             'NKR': 'X',
             'NKS': 'X',
             'NKT': 'X',
             'NKV': 'X',
             'NKW': 'X',
             'NKY': 'X',
             'NMA': 'X',
             'NMB': 'X',
             'NMC': 'X',
             'NMD': 'X',
             'NMG': 'X',
             'NMH': 'X',
             'NMK': 'X',
             'NMM': 'X',
             'NMN': 'X',
             'NMR': 'X',
             'NMS': 'X',
             'NMT': 'X',
             'NMV': 'X',
             'NMW': 'X',
             'NMY': 'X',
             'NNA': 'X',
             'NNB': 'X',
             'NNC': 'X',
             'NND': 'X',
             'NNG': 'X',
             'NNH': 'X',
             'NNK': 'X',
             'NNM': 'X',
             'NNN': 'X',
             'NNR': 'X',
             'NNS': 'X',
             'NNT': 'X',
             'NNV': 'X',
             'NNW': 'X',
             'NNY': 'X',
             'NRA': 'X',
             'NRB': 'X',
             'NRC': 'X',
             'NRD': 'X',
             'NRG': 'X',
             'NRH': 'X',
             'NRK': 'X',
             'NRM': 'X',
             'NRN': 'X',
             'NRR': 'X',
             'NRS': 'X',
             'NRT': 'X',
             'NRV': 'X',
             'NRW': 'X',
             'NRY': 'X',
             'NSA': 'X',
             'NSB': 'X',
             'NSC': 'X',
             'NSD': 'X',
             'NSG': 'X',
             'NSH': 'X',
             'NSK': 'X',
             'NSM': 'X',
             'NSN': 'X',
             'NSR': 'X',
             'NSS': 'X',
             'NST': 'X',
             'NSV': 'X',
             'NSW': 'X',
             'NSY': 'X',
             'NTA': 'ILV',
             'NTB': 'X',
             'NTC': 'ILVF',
             'NTD': 'X',
             'NTG': 'MLV',
             'NTH': 'ILVF',
             'NTK': 'X',
             'NTM': 'ILVF',
             'NTN': 'X',
             'NTR': 'IMLV',
             'NTS': 'X',
             'NTT': 'ILVF',
             'NTV': 'X',
             'NTW': 'ILVF',
             'NTY': 'ILVF',
             'NVA': 'X',
             'NVB': 'X',
             'NVC': 'X',
             'NVD': 'X',
             'NVG': 'X',
             'NVH': 'X',
             'NVK': 'X',
             'NVM': 'X',
             'NVN': 'X',
             'NVR': 'X',
             'NVS': 'X',
             'NVT': 'X',
             'NVV': 'X',
             'NVW': 'X',
             'NVY': 'X',
             'NWA': 'X',
             'NWB': 'X',
             'NWC': 'X',
             'NWD': 'X',
             'NWG': 'X',
             'NWH': 'X',
             'NWK': 'X',
             'NWM': 'X',
             'NWN': 'X',
             'NWR': 'X',
             'NWS': 'X',
             'NWT': 'X',
             'NWV': 'X',
             'NWW': 'X',
             'NWY': 'X',
             'NYA': 'X',
             'NYB': 'X',
             'NYC': 'X',
             'NYD': 'X',
             'NYG': 'X',
             'NYH': 'X',
             'NYK': 'X',
             'NYM': 'X',
             'NYN': 'X',
             'NYR': 'X',
             'NYS': 'X',
             'NYT': 'X',
             'NYV': 'X',
             'NYW': 'X',
             'NYY': 'X',
             'RAA': 'KE',
             'RAB': 'NKDE',
             'RAC': 'ND',
             'RAD': 'KNED',
             'RAG': 'KE',
             'RAH': 'KNED',
             'RAK': 'KNED',
             'RAM': 'KNED',
             'RAN': 'KNED',
             'RAR': 'KE',
             'RAS': 'NKDE',
             'RAT': 'ND',
             'RAV': 'KNED',
             'RAW': 'KNED',
             'RAY': 'ND',
             'RBA': 'X',
             'RBB': 'X',
             'RBC': 'X',
             'RBD': 'X',
             'RBG': 'X',
             'RBH': 'X',
             'RBK': 'X',
             'RBM': 'X',
             'RBN': 'X',
             'RBR': 'X',
             'RBS': 'X',
             'RBT': 'X',
             'RBV': 'X',
             'RBW': 'X',
             'RBY': 'X',
             'RCA': 'TA',
             'RCB': 'TA',
             'RCC': 'TA',
             'RCD': 'TA',
             'RCG': 'TA',
             'RCH': 'TA',
             'RCK': 'TA',
             'RCM': 'TA',
             'RCN': 'TA',
             'RCR': 'TA',
             'RCS': 'TA',
             'RCT': 'TA',
             'RCV': 'TA',
             'RCW': 'TA',
             'RCY': 'TA',
             'RDA': 'X',
             'RDB': 'X',
             'RDC': 'X',
             'RDD': 'X',
             'RDG': 'X',
             'RDH': 'X',
             'RDK': 'X',
             'RDM': 'X',
             'RDN': 'X',
             'RDR': 'X',
             'RDS': 'X',
             'RDT': 'X',
             'RDV': 'X',
             'RDW': 'X',
             'RDY': 'X',
             'RGA': 'RG',
             'RGB': 'SRG',
             'RGC': 'SG',
             'RGD': 'RSG',
             'RGG': 'RG',
             'RGH': 'RSG',
             'RGK': 'RSG',
             'RGM': 'RSG',
             'RGN': 'RSG',
             'RGR': 'RG',
             'RGS': 'SRG',
             'RGT': 'SG',
             'RGV': 'RSG',
             'RGW': 'RSG',
             'RGY': 'SG',
             'RHA': 'X',
             'RHB': 'X',
             'RHC': 'X',
             'RHD': 'X',
             'RHG': 'X',
             'RHH': 'X',
             'RHK': 'X',
             'RHM': 'X',
             'RHN': 'X',
             'RHR': 'X',
             'RHS': 'X',
             'RHT': 'X',
             'RHV': 'X',
             'RHW': 'X',
             'RHY': 'X',
             'RKA': 'RIGV',
             'RKB': 'X',
             'RKC': 'SIGV',
             'RKD': 'X',
             'RKG': 'RMGV',
             'RKH': 'X',
             'RKK': 'X',
             'RKM': 'X',
             'RKN': 'X',
             'RKR': 'X',
             'RKS': 'X',
             'RKT': 'SIGV',
             'RKV': 'X',
             'RKW': 'X',
             'RKY': 'SIGV',
             'RMA': 'KTEA',
             'RMB': 'X',
             'RMC': 'NTDA',
             'RMD': 'X',
             'RMG': 'KTEA',
             'RMH': 'X',
             'RMK': 'X',
             'RMM': 'X',
             'RMN': 'X',
             'RMR': 'KTEA',
             'RMS': 'X',
             'RMT': 'NTDA',
             'RMV': 'X',
             'RMW': 'X',
             'RMY': 'NTDA',
             'RNA': 'X',
             'RNB': 'X',
             'RNC': 'X',
             'RND': 'X',
             'RNG': 'X',
             'RNH': 'X',
             'RNK': 'X',
             'RNM': 'X',
             'RNN': 'X',
             'RNR': 'X',
             'RNS': 'X',
             'RNT': 'X',
             'RNV': 'X',
             'RNW': 'X',
             'RNY': 'X',
             'RRA': 'KREG',
             'RRB': 'X',
             'RRC': 'NSDG',
             'RRD': 'X',
             'RRG': 'KREG',
             'RRH': 'X',
             'RRK': 'X',
             'RRM': 'X',
             'RRN': 'X',
             'RRR': 'KREG',
             'RRS': 'X',
             'RRT': 'NSDG',
             'RRV': 'X',
             'RRW': 'X',
             'RRY': 'NSDG',
             'RSA': 'TRAG',
             'RSB': 'X',
             'RSC': 'TSAG',
             'RSD': 'X',
             'RSG': 'TRAG',
             'RSH': 'X',
             'RSK': 'X',
             'RSM': 'X',
             'RSN': 'X',
             'RSR': 'TRAG',
             'RSS': 'X',
             'RST': 'TSAG',
             'RSV': 'X',
             'RSW': 'X',
             'RSY': 'TSAG',
             'RTA': 'IV',
             'RTB': 'IMV',
             'RTC': 'IV',
             'RTD': 'IMV',
             'RTG': 'MV',
             'RTH': 'IV',
             'RTK': 'MIV',
             'RTM': 'IV',
             'RTN': 'IMV',
             'RTR': 'IMV',
             'RTS': 'IMV',
             'RTT': 'IV',
             'RTV': 'IMV',
             'RTW': 'IV',
             'RTY': 'IV',
             'RVA': 'X',
             'RVB': 'X',
             'RVC': 'X',
             'RVD': 'X',
             'RVG': 'X',
             'RVH': 'X',
             'RVK': 'X',
             'RVM': 'X',
             'RVN': 'X',
             'RVR': 'X',
             'RVS': 'X',
             'RVT': 'X',
             'RVV': 'X',
             'RVW': 'X',
             'RVY': 'X',
             'RWA': 'KIEV',
             'RWB': 'X',
             'RWC': 'NIDV',
             'RWD': 'X',
             'RWG': 'KMEV',
             'RWH': 'X',
             'RWK': 'X',
             'RWM': 'X',
             'RWN': 'X',
             'RWR': 'X',
             'RWS': 'X',
             'RWT': 'NIDV',
             'RWV': 'X',
             'RWW': 'X',
             'RWY': 'NIDV',
             'RYA': 'TIAV',
             'RYB': 'X',
             'RYC': 'TIAV',
             'RYD': 'X',
             'RYG': 'TMAV',
             'RYH': 'TIAV',
             'RYK': 'X',
             'RYM': 'TIAV',
             'RYN': 'X',
             'RYR': 'X',
             'RYS': 'X',
             'RYT': 'TIAV',
             'RYV': 'X',
             'RYW': 'TIAV',
             'RYY': 'TIAV',
             'SAA': 'QE',
             'SAB': 'HQDE',
             'SAC': 'HD',
             'SAD': 'QHED',
             'SAG': 'QE',
             'SAH': 'QHED',
             'SAK': 'QHED',
             'SAM': 'QHED',
             'SAN': 'QHED',
             'SAR': 'QE',
             'SAS': 'HQDE',
             'SAT': 'HD',
             'SAV': 'QHED',
             'SAW': 'QHED',
             'SAY': 'HD',
             'SBA': 'X',
             'SBB': 'X',
             'SBC': 'X',
             'SBD': 'X',
             'SBG': 'X',
             'SBH': 'X',
             'SBK': 'X',
             'SBM': 'X',
             'SBN': 'X',
             'SBR': 'X',
             'SBS': 'X',
             'SBT': 'X',
             'SBV': 'X',
             'SBW': 'X',
             'SBY': 'X',
             'SCA': 'PA',
             'SCB': 'PA',
             'SCC': 'PA',
             'SCD': 'PA',
             'SCG': 'PA',
             'SCH': 'PA',
             'SCK': 'PA',
             'SCM': 'PA',
             'SCN': 'PA',
             'SCR': 'PA',
             'SCS': 'PA',
             'SCT': 'PA',
             'SCV': 'PA',
             'SCW': 'PA',
             'SCY': 'PA',
             'SDA': 'X',
             'SDB': 'X',
             'SDC': 'X',
             'SDD': 'X',
             'SDG': 'X',
             'SDH': 'X',
             'SDK': 'X',
             'SDM': 'X',
             'SDN': 'X',
             'SDR': 'X',
             'SDS': 'X',
             'SDT': 'X',
             'SDV': 'X',
             'SDW': 'X',
             'SDY': 'X',
             'SGA': 'RG',
             'SGB': 'RG',
             'SGC': 'RG',
             'SGD': 'RG',
             'SGG': 'RG',
             'SGH': 'RG',
             'SGK': 'RG',
             'SGM': 'RG',
             'SGN': 'RG',
             'SGR': 'RG',
             'SGS': 'RG',
             'SGT': 'RG',
             'SGV': 'RG',
             'SGW': 'RG',
             'SGY': 'RG',
             'SHA': 'X',
             'SHB': 'X',
             'SHC': 'X',
             'SHD': 'X',
             'SHG': 'X',
             'SHH': 'X',
             'SHK': 'X',
             'SHM': 'X',
             'SHN': 'X',
             'SHR': 'X',
             'SHS': 'X',
             'SHT': 'X',
             'SHV': 'X',
             'SHW': 'X',
             'SHY': 'X',
             'SKA': 'RLGV',
             'SKB': 'RLGV',
             'SKC': 'RLGV',
             'SKD': 'RLGV',
             'SKG': 'RLGV',
             'SKH': 'RLGV',
             'SKK': 'RLGV',
             'SKM': 'RLGV',
             'SKN': 'RLGV',
             'SKR': 'RLGV',
             'SKS': 'RLGV',
             'SKT': 'RLGV',
             'SKV': 'RLGV',
             'SKW': 'RLGV',
             'SKY': 'RLGV',
             'SMA': 'QPEA',
             'SMB': 'X',
             'SMC': 'HPDA',
             'SMD': 'X',
             'SMG': 'QPEA',
             'SMH': 'X',
             'SMK': 'X',
             'SMM': 'X',
             'SMN': 'X',
             'SMR': 'QPEA',
             'SMS': 'X',
             'SMT': 'HPDA',
             'SMV': 'X',
             'SMW': 'X',
             'SMY': 'HPDA',
             'SNA': 'X',
             'SNB': 'X',
             'SNC': 'X',
             'SND': 'X',
             'SNG': 'X',
             'SNH': 'X',
             'SNK': 'X',
             'SNM': 'X',
             'SNN': 'X',
             'SNR': 'X',
             'SNS': 'X',
             'SNT': 'X',
             'SNV': 'X',
             'SNW': 'X',
             'SNY': 'X',
             'SRA': 'QREG',
             'SRB': 'X',
             'SRC': 'HRDG',
             'SRD': 'X',
             'SRG': 'QREG',
             'SRH': 'X',
             'SRK': 'X',
             'SRM': 'X',
             'SRN': 'X',
             'SRR': 'QREG',
             'SRS': 'X',
             'SRT': 'HRDG',
             'SRV': 'X',
             'SRW': 'X',
             'SRY': 'HRDG',
             'SSA': 'PRAG',
             'SSB': 'PRAG',
             'SSC': 'PRAG',
             'SSD': 'PRAG',
             'SSG': 'PRAG',
             'SSH': 'PRAG',
             'SSK': 'PRAG',
             'SSM': 'PRAG',
             'SSN': 'PRAG',
             'SSR': 'PRAG',
             'SSS': 'PRAG',
             'SST': 'PRAG',
             'SSV': 'PRAG',
             'SSW': 'PRAG',
             'SSY': 'PRAG',
             'STA': 'LV',
             'STB': 'LV',
             'STC': 'LV',
             'STD': 'LV',
             'STG': 'LV',
             'STH': 'LV',
             'STK': 'LV',
             'STM': 'LV',
             'STN': 'LV',
             'STR': 'LV',
             'STS': 'LV',
             'STT': 'LV',
             'STV': 'LV',
             'STW': 'LV',
             'STY': 'LV',
             'SVA': 'X',
             'SVB': 'X',
             'SVC': 'X',
             'SVD': 'X',
             'SVG': 'X',
             'SVH': 'X',
             'SVK': 'X',
             'SVM': 'X',
             'SVN': 'X',
             'SVR': 'X',
             'SVS': 'X',
             'SVT': 'X',
             'SVV': 'X',
             'SVW': 'X',
             'SVY': 'X',
             'SWA': 'QLEV',
             'SWB': 'X',
             'SWC': 'HLDV',
             'SWD': 'X',
             'SWG': 'QLEV',
             'SWH': 'X',
             'SWK': 'X',
             'SWM': 'X',
             'SWN': 'X',
             'SWR': 'QLEV',
             'SWS': 'X',
             'SWT': 'HLDV',
             'SWV': 'X',
             'SWW': 'X',
             'SWY': 'HLDV',
             'SYA': 'PLAV',
             'SYB': 'PLAV',
             'SYC': 'PLAV',
             'SYD': 'PLAV',
             'SYG': 'PLAV',
             'SYH': 'PLAV',
             'SYK': 'PLAV',
             'SYM': 'PLAV',
             'SYN': 'PLAV',
             'SYR': 'PLAV',
             'SYS': 'PLAV',
             'SYT': 'PLAV',
             'SYV': 'PLAV',
             'SYW': 'PLAV',
             'SYY': 'PLAV',
             'TAA': '*',
             'TAB': 'Y*',
             'TAC': 'Y',
             'TAD': '*Y',
             'TAG': '*',
             'TAH': '*Y',
             'TAK': '*Y',
             'TAM': '*Y',
             'TAN': '*Y',
             'TAR': '*',
             'TAS': 'Y*',
             'TAT': 'Y',
             'TAV': '*Y',
             'TAW': '*Y',
             'TAY': 'Y',
             'TBA': 'S*L',
             'TBB': 'X',
             'TBC': 'SCF',
             'TBD': 'X',
             'TBG': 'SWL',
             'TBH': 'X',
             'TBK': 'X',
             'TBM': 'X',
             'TBN': 'X',
             'TBR': 'S*WL',
             'TBS': 'X',
             'TBT': 'SCF',
             'TBV': 'X',
             'TBW': 'X',
             'TBY': 'SCF',
             'TCA': 'S',
             'TCB': 'S',
             'TCC': 'S',
             'TCD': 'S',
             'TCG': 'S',
             'TCH': 'S',
             'TCK': 'S',
             'TCM': 'S',
             'TCN': 'S',
             'TCR': 'S',
             'TCS': 'S',
             'TCT': 'S',
             'TCV': 'S',
             'TCW': 'S',
             'TCY': 'S',
             'TDA': '*L',
             'TDB': 'X',
             'TDC': 'YCF',
             'TDD': 'X',
             'TDG': '*WL',
             'TDH': 'X',
             'TDK': 'X',
             'TDM': 'X',
             'TDN': 'X',
             'TDR': '*WL',
             'TDS': 'X',
             'TDT': 'YCF',
             'TDV': 'X',
             'TDW': 'X',
             'TDY': 'YCF',
             'TGA': '*',
             'TGB': 'CW',
             'TGC': 'C',
             'TGD': '*WC',
             'TGG': 'W',
             'TGH': '*C',
             'TGK': 'WC',
             'TGM': '*C',
             'TGN': '*CW',
             'TGR': '*W',
             'TGS': 'CW',
             'TGT': 'C',
             'TGV': '*CW',
             'TGW': '*C',
             'TGY': 'C',
             'THA': '*SL',
             'THB': 'X',
             'THC': 'YSF',
             'THD': 'X',
             'THG': '*SL',
             'THH': 'X',
             'THK': 'X',
             'THM': 'X',
             'THN': 'X',
             'THR': '*SL',
             'THS': 'X',
             'THT': 'YSF',
             'THV': 'X',
             'THW': 'X',
             'THY': 'YSF',
             'TKA': '*L',
             'TKB': 'CWFL',
             'TKC': 'CF',
             'TKD': 'X',
             'TKG': 'WL',
             'TKH': '*CLF',
             'TKK': 'WCLF',
             'TKM': '*CLF',
             'TKN': 'X',
             'TKR': '*WL',
             'TKS': 'CWFL',
             'TKT': 'CF',
             'TKV': 'X',
             'TKW': '*CLF',
             'TKY': 'CF',
             'TMA': '*S',
             'TMB': 'Y*S',
             'TMC': 'YS',
             'TMD': '*YS',
             'TMG': '*S',
             'TMH': '*YS',
             'TMK': '*YS',
             'TMM': '*YS',
             'TMN': '*YS',
             'TMR': '*S',
             'TMS': 'Y*S',
             'TMT': 'YS',
             'TMV': '*YS',
             'TMW': '*YS',
             'TMY': 'YS',
             'TNA': '*SL',
             'TNB': 'X',
             'TNC': 'YSCF',
             'TND': 'X',
             'TNG': '*SWL',
             'TNH': 'X',
             'TNK': 'X',
             'TNM': 'X',
             'TNN': 'X',
             'TNR': '*SWL',
             'TNS': 'X',
             'TNT': 'YSCF',
             'TNV': 'X',
             'TNW': 'X',
             'TNY': 'YSCF',
             'TRA': '*',
             'TRB': 'Y*CW',
             'TRC': 'YC',
             'TRD': '*YWC',
             'TRG': '*W',
             'TRH': '*YC',
             'TRK': '*YWC',
             'TRM': '*YC',
             'TRN': '*YCW',
             'TRR': '*W',
             'TRS': 'Y*CW',
             'TRT': 'YC',
             'TRV': '*YCW',
             'TRW': '*YC',
             'TRY': 'YC',
             'TSA': 'S*',
             'TSB': 'SCW',
             'TSC': 'SC',
             'TSD': 'S*WC',
             'TSG': 'SW',
             'TSH': 'S*C',
             'TSK': 'SWC',
             'TSM': 'S*C',
             'TSN': 'S*CW',
             'TSR': 'S*W',
             'TSS': 'SCW',
             'TST': 'SC',
             'TSV': 'S*CW',
             'TSW': 'S*C',
             'TSY': 'SC',
             'TTA': 'L',
             'TTB': 'FL',
             'TTC': 'F',
             'TTD': 'LF',
             'TTG': 'L',
             'TTH': 'LF',
             'TTK': 'LF',
             'TTM': 'LF',
             'TTN': 'LF',
             'TTR': 'L',
             'TTS': 'FL',
             'TTT': 'F',
             'TTV': 'LF',
             'TTW': 'LF',
             'TTY': 'F',
             'TVA': '*S',
             'TVB': 'X',
             'TVC': 'YSC',
             'TVD': 'X',
             'TVG': '*SW',
             'TVH': '*YSC',
             'TVK': 'X',
             'TVM': '*YSC',
             'TVN': 'X',
             'TVR': '*SW',
             'TVS': 'X',
             'TVT': 'YSC',
             'TVV': 'X',
             'TVW': '*YSC',
             'TVY': 'YSC',
             'TWA': '*L',
             'TWB': 'Y*FL',
             'TWC': 'YF',
             'TWD': '*YLF',
             'TWG': '*L',
             'TWH': '*YLF',
             'TWK': '*YLF',
             'TWM': '*YLF',
             'TWN': '*YLF',
             'TWR': '*L',
             'TWS': 'Y*FL',
             'TWT': 'YF',
             'TWV': '*YLF',
             'TWW': '*YLF',
             'TWY': 'YF',
             'TYA': 'SL',
             'TYB': 'SFL',
             'TYC': 'SF',
             'TYD': 'SLF',
             'TYG': 'SL',
             'TYH': 'SLF',
             'TYK': 'SLF',
             'TYM': 'SLF',
             'TYN': 'SLF',
             'TYR': 'SL',
             'TYS': 'SFL',
             'TYT': 'SF',
             'TYV': 'SLF',
             'TYW': 'SLF',
             'TYY': 'SF',
             'VAA': 'KQE',
             'VAB': 'X',
             'VAC': 'NHD',
             'VAD': 'X',
             'VAG': 'KQE',
             'VAH': 'X',
             'VAK': 'X',
             'VAM': 'X',
             'VAN': 'X',
             'VAR': 'KQE',
             'VAS': 'X',
             'VAT': 'NHD',
             'VAV': 'X',
             'VAW': 'X',
             'VAY': 'NHD',
             'VBA': 'X',
             'VBB': 'X',
             'VBC': 'X',
             'VBD': 'X',
             'VBG': 'X',
             'VBH': 'X',
             'VBK': 'X',
             'VBM': 'X',
             'VBN': 'X',
             'VBR': 'X',
             'VBS': 'X',
             'VBT': 'X',
             'VBV': 'X',
             'VBW': 'X',
             'VBY': 'X',
             'VCA': 'TPA',
             'VCB': 'TPA',
             'VCC': 'TPA',
             'VCD': 'TPA',
             'VCG': 'TPA',
             'VCH': 'TPA',
             'VCK': 'TPA',
             'VCM': 'TPA',
             'VCN': 'TPA',
             'VCR': 'TPA',
             'VCS': 'TPA',
             'VCT': 'TPA',
             'VCV': 'TPA',
             'VCW': 'TPA',
             'VCY': 'TPA',
             'VDA': 'X',
             'VDB': 'X',
             'VDC': 'X',
             'VDD': 'X',
             'VDG': 'X',
             'VDH': 'X',
             'VDK': 'X',
             'VDM': 'X',
             'VDN': 'X',
             'VDR': 'X',
             'VDS': 'X',
             'VDT': 'X',
             'VDV': 'X',
             'VDW': 'X',
             'VDY': 'X',
             'VGA': 'RG',
             'VGB': 'SRG',
             'VGC': 'SRG',
             'VGD': 'RSG',
             'VGG': 'RG',
             'VGH': 'RSG',
             'VGK': 'RSG',
             'VGM': 'RSG',
             'VGN': 'RSG',
             'VGR': 'RG',
             'VGS': 'SRG',
             'VGT': 'SRG',
             'VGV': 'RSG',
             'VGW': 'RSG',
             'VGY': 'SRG',
             'VHA': 'X',
             'VHB': 'X',
             'VHC': 'X',
             'VHD': 'X',
             'VHG': 'X',
             'VHH': 'X',
             'VHK': 'X',
             'VHM': 'X',
             'VHN': 'X',
             'VHR': 'X',
             'VHS': 'X',
             'VHT': 'X',
             'VHV': 'X',
             'VHW': 'X',
             'VHY': 'X',
             'VKA': 'X',
             'VKB': 'X',
             'VKC': 'X',
             'VKD': 'X',
             'VKG': 'X',
             'VKH': 'X',
             'VKK': 'X',
             'VKM': 'X',
             'VKN': 'X',
             'VKR': 'X',
             'VKS': 'X',
             'VKT': 'X',
             'VKV': 'X',
             'VKW': 'X',
             'VKY': 'X',
             'VMA': 'X',
             'VMB': 'X',
             'VMC': 'X',
             'VMD': 'X',
             'VMG': 'X',
             'VMH': 'X',
             'VMK': 'X',
             'VMM': 'X',
             'VMN': 'X',
             'VMR': 'X',
             'VMS': 'X',
             'VMT': 'X',
             'VMV': 'X',
             'VMW': 'X',
             'VMY': 'X',
             'VNA': 'X',
             'VNB': 'X',
             'VNC': 'X',
             'VND': 'X',
             'VNG': 'X',
             'VNH': 'X',
             'VNK': 'X',
             'VNM': 'X',
             'VNN': 'X',
             'VNR': 'X',
             'VNS': 'X',
             'VNT': 'X',
             'VNV': 'X',
             'VNW': 'X',
             'VNY': 'X',
             'VRA': 'X',
             'VRB': 'X',
             'VRC': 'X',
             'VRD': 'X',
             'VRG': 'X',
             'VRH': 'X',
             'VRK': 'X',
             'VRM': 'X',
             'VRN': 'X',
             'VRR': 'X',
             'VRS': 'X',
             'VRT': 'X',
             'VRV': 'X',
             'VRW': 'X',
             'VRY': 'X',
             'VSA': 'X',
             'VSB': 'X',
             'VSC': 'X',
             'VSD': 'X',
             'VSG': 'X',
             'VSH': 'X',
             'VSK': 'X',
             'VSM': 'X',
             'VSN': 'X',
             'VSR': 'X',
             'VSS': 'X',
             'VST': 'X',
             'VSV': 'X',
             'VSW': 'X',
             'VSY': 'X',
             'VTA': 'ILV',
             'VTB': 'IMLV',
             'VTC': 'ILV',
             'VTD': 'IMLV',
             'VTG': 'MLV',
             'VTH': 'ILV',
             'VTK': 'MILV',
             'VTM': 'ILV',
             'VTN': 'IMLV',
             'VTR': 'IMLV',
             'VTS': 'IMLV',
             'VTT': 'ILV',
             'VTV': 'IMLV',
             'VTW': 'ILV',
             'VTY': 'ILV',
             'VVA': 'X',
             'VVB': 'X',
             'VVC': 'X',
             'VVD': 'X',
             'VVG': 'X',
             'VVH': 'X',
             'VVK': 'X',
             'VVM': 'X',
             'VVN': 'X',
             'VVR': 'X',
             'VVS': 'X',
             'VVT': 'X',
             'VVV': 'X',
             'VVW': 'X',
             'VVY': 'X',
             'VWA': 'X',
             'VWB': 'X',
             'VWC': 'X',
             'VWD': 'X',
             'VWG': 'X',
             'VWH': 'X',
             'VWK': 'X',
             'VWM': 'X',
             'VWN': 'X',
             'VWR': 'X',
             'VWS': 'X',
             'VWT': 'X',
             'VWV': 'X',
             'VWW': 'X',
             'VWY': 'X',
             'VYA': 'X',
             'VYB': 'X',
             'VYC': 'X',
             'VYD': 'X',
             'VYG': 'X',
             'VYH': 'X',
             'VYK': 'X',
             'VYM': 'X',
             'VYN': 'X',
             'VYR': 'X',
             'VYS': 'X',
             'VYT': 'X',
             'VYV': 'X',
             'VYW': 'X',
             'VYY': 'X',
             'WAA': 'K*',
             'WAB': 'NKY*',
             'WAC': 'NY',
             'WAD': 'KN*Y',
             'WAG': 'K*',
             'WAH': 'KN*Y',
             'WAK': 'KN*Y',
             'WAM': 'KN*Y',
             'WAN': 'KN*Y',
             'WAR': 'K*',
             'WAS': 'NKY*',
             'WAT': 'NY',
             'WAV': 'KN*Y',
             'WAW': 'KN*Y',
             'WAY': 'NY',
             'WBA': 'X',
             'WBB': 'X',
             'WBC': 'X',
             'WBD': 'X',
             'WBG': 'X',
             'WBH': 'X',
             'WBK': 'X',
             'WBM': 'X',
             'WBN': 'X',
             'WBR': 'X',
             'WBS': 'X',
             'WBT': 'X',
             'WBV': 'X',
             'WBW': 'X',
             'WBY': 'X',
             'WCA': 'TS',
             'WCB': 'TS',
             'WCC': 'TS',
             'WCD': 'TS',
             'WCG': 'TS',
             'WCH': 'TS',
             'WCK': 'TS',
             'WCM': 'TS',
             'WCN': 'TS',
             'WCR': 'TS',
             'WCS': 'TS',
             'WCT': 'TS',
             'WCV': 'TS',
             'WCW': 'TS',
             'WCY': 'TS',
             'WDA': 'X',
             'WDB': 'X',
             'WDC': 'X',
             'WDD': 'X',
             'WDG': 'X',
             'WDH': 'X',
             'WDK': 'X',
             'WDM': 'X',
             'WDN': 'X',
             'WDR': 'X',
             'WDS': 'X',
             'WDT': 'X',
             'WDV': 'X',
             'WDW': 'X',
             'WDY': 'X',
             'WGA': 'R*',
             'WGB': 'SRCW',
             'WGC': 'SC',
             'WGD': 'X',
             'WGG': 'RW',
             'WGH': 'RS*C',
             'WGK': 'RSWC',
             'WGM': 'RS*C',
             'WGN': 'X',
             'WGR': 'R*W',
             'WGS': 'SRCW',
             'WGT': 'SC',
             'WGV': 'X',
             'WGW': 'RS*C',
             'WGY': 'SC',
             'WHA': 'X',
             'WHB': 'X',
             'WHC': 'X',
             'WHD': 'X',
             'WHG': 'X',
             'WHH': 'X',
             'WHK': 'X',
             'WHM': 'X',
             'WHN': 'X',
             'WHR': 'X',
             'WHS': 'X',
             'WHT': 'X',
             'WHV': 'X',
             'WHW': 'X',
             'WHY': 'X',
             'WKA': 'RI*L',
             'WKB': 'X',
             'WKC': 'SICF',
             'WKD': 'X',
             'WKG': 'RMWL',
             'WKH': 'X',
             'WKK': 'X',
             'WKM': 'X',
             'WKN': 'X',
             'WKR': 'X',
             'WKS': 'X',
             'WKT': 'SICF',
             'WKV': 'X',
             'WKW': 'X',
             'WKY': 'SICF',
             'WMA': 'KT*S',
             'WMB': 'X',
             'WMC': 'NTYS',
             'WMD': 'X',
             'WMG': 'KT*S',
             'WMH': 'X',
             'WMK': 'X',
             'WMM': 'X',
             'WMN': 'X',
             'WMR': 'KT*S',
             'WMS': 'X',
             'WMT': 'NTYS',
             'WMV': 'X',
             'WMW': 'X',
             'WMY': 'NTYS',
             'WNA': 'X',
             'WNB': 'X',
             'WNC': 'X',
             'WND': 'X',
             'WNG': 'X',
             'WNH': 'X',
             'WNK': 'X',
             'WNM': 'X',
             'WNN': 'X',
             'WNR': 'X',
             'WNS': 'X',
             'WNT': 'X',
             'WNV': 'X',
             'WNW': 'X',
             'WNY': 'X',
             'WRA': 'KR*',
             'WRB': 'X',
             'WRC': 'NSYC',
             'WRD': 'X',
             'WRG': 'KR*W',
             'WRH': 'X',
             'WRK': 'X',
             'WRM': 'X',
             'WRN': 'X',
             'WRR': 'KR*W',
             'WRS': 'X',
             'WRT': 'NSYC',
             'WRV': 'X',
             'WRW': 'X',
             'WRY': 'NSYC',
             'WSA': 'TRS*',
             'WSB': 'X',
             'WSC': 'TSC',
             'WSD': 'X',
             'WSG': 'TRSW',
             'WSH': 'X',
             'WSK': 'X',
             'WSM': 'X',
             'WSN': 'X',
             'WSR': 'X',
             'WSS': 'X',
             'WST': 'TSC',
             'WSV': 'X',
             'WSW': 'X',
             'WSY': 'TSC',
             'WTA': 'IL',
             'WTB': 'IMFL',
             'WTC': 'IF',
             'WTD': 'IMLF',
             'WTG': 'ML',
             'WTH': 'ILF',
             'WTK': 'MILF',
             'WTM': 'ILF',
             'WTN': 'IMLF',
             'WTR': 'IML',
             'WTS': 'IMFL',
             'WTT': 'IF',
             'WTV': 'IMLF',
             'WTW': 'ILF',
             'WTY': 'IF',
             'WVA': 'X',
             'WVB': 'X',
             'WVC': 'X',
             'WVD': 'X',
             'WVG': 'X',
             'WVH': 'X',
             'WVK': 'X',
             'WVM': 'X',
             'WVN': 'X',
             'WVR': 'X',
             'WVS': 'X',
             'WVT': 'X',
             'WVV': 'X',
             'WVW': 'X',
             'WVY': 'X',
             'WWA': 'KI*L',
             'WWB': 'X',
             'WWC': 'NIYF',
             'WWD': 'X',
             'WWG': 'KM*L',
             'WWH': 'X',
             'WWK': 'X',
             'WWM': 'X',
             'WWN': 'X',
             'WWR': 'X',
             'WWS': 'X',
             'WWT': 'NIYF',
             'WWV': 'X',
             'WWW': 'X',
             'WWY': 'NIYF',
             'WYA': 'TISL',
             'WYB': 'X',
             'WYC': 'TISF',
             'WYD': 'X',
             'WYG': 'TMSL',
             'WYH': 'X',
             'WYK': 'X',
             'WYM': 'X',
             'WYN': 'X',
             'WYR': 'X',
             'WYS': 'X',
             'WYT': 'TISF',
             'WYV': 'X',
             'WYW': 'X',
             'WYY': 'TISF',
             'YAA': 'Q*',
             'YAB': 'HQY*',
             'YAC': 'HY',
             'YAD': 'QH*Y',
             'YAG': 'Q*',
             'YAH': 'QH*Y',
             'YAK': 'QH*Y',
             'YAM': 'QH*Y',
             'YAN': 'QH*Y',
             'YAR': 'Q*',
             'YAS': 'HQY*',
             'YAT': 'HY',
             'YAV': 'QH*Y',
             'YAW': 'QH*Y',
             'YAY': 'HY',
             'YBA': 'X',
             'YBB': 'X',
             'YBC': 'X',
             'YBD': 'X',
             'YBG': 'X',
             'YBH': 'X',
             'YBK': 'X',
             'YBM': 'X',
             'YBN': 'X',
             'YBR': 'X',
             'YBS': 'X',
             'YBT': 'X',
             'YBV': 'X',
             'YBW': 'X',
             'YBY': 'X',
             'YCA': 'PS',
             'YCB': 'PS',
             'YCC': 'PS',
             'YCD': 'PS',
             'YCG': 'PS',
             'YCH': 'PS',
             'YCK': 'PS',
             'YCM': 'PS',
             'YCN': 'PS',
             'YCR': 'PS',
             'YCS': 'PS',
             'YCT': 'PS',
             'YCV': 'PS',
             'YCW': 'PS',
             'YCY': 'PS',
             'YDA': 'QRL*',
             'YDB': 'X',
             'YDC': 'X',
             'YDD': 'X',
             'YDG': 'X',
             'YDH': 'X',
             'YDK': 'X',
             'YDM': 'X',
             'YDN': 'X',
             'YDR': 'X',
             'YDS': 'X',
             'YDT': 'X',
             'YDV': 'X',
             'YDW': 'X',
             'YDY': 'X',
             'YGA': 'R*',
             'YGB': 'RCW',
             'YGC': 'RC',
             'YGD': 'R*WC',
             'YGG': 'RW',
             'YGH': 'R*C',
             'YGK': 'RWC',
             'YGM': 'R*C',
             'YGN': 'R*CW',
             'YGR': 'R*W',
             'YGS': 'RCW',
             'YGT': 'RC',
             'YGV': 'R*CW',
             'YGW': 'R*C',
             'YGY': 'RC',
             'YHA': 'X',
             'YHB': 'X',
             'YHC': 'X',
             'YHD': 'X',
             'YHG': 'X',
             'YHH': 'X',
             'YHK': 'X',
             'YHM': 'X',
             'YHN': 'X',
             'YHR': 'X',
             'YHS': 'X',
             'YHT': 'X',
             'YHV': 'X',
             'YHW': 'X',
             'YHY': 'X',
             'YKA': 'RL*',
             'YKB': 'X',
             'YKC': 'RLCF',
             'YKD': 'X',
             'YKG': 'RLW',
             'YKH': 'X',
             'YKK': 'X',
             'YKM': 'X',
             'YKN': 'X',
             'YKR': 'RL*W',
             'YKS': 'X',
             'YKT': 'RLCF',
             'YKV': 'X',
             'YKW': 'X',
             'YKY': 'RLCF',
             'YMA': 'QP*S',
             'YMB': 'X',
             'YMC': 'HPYS',
             'YMD': 'X',
             'YMG': 'QP*S',
             'YMH': 'X',
             'YMK': 'X',
             'YMM': 'X',
             'YMN': 'X',
             'YMR': 'QP*S',
             'YMS': 'X',
             'YMT': 'HPYS',
             'YMV': 'X',
             'YMW': 'X',
             'YMY': 'HPYS',
             'YNA': 'X',
             'YNB': 'X',
             'YNC': 'X',
             'YND': 'X',
             'YNG': 'X',
             'YNH': 'X',
             'YNK': 'X',
             'YNM': 'X',
             'YNN': 'X',
             'YNR': 'X',
             'YNS': 'X',
             'YNT': 'X',
             'YNV': 'X',
             'YNW': 'X',
             'YNY': 'X',
             'YRA': 'QR*',
             'YRB': 'X',
             'YRC': 'HRYC',
             'YRD': 'X',
             'YRG': 'QR*W',
             'YRH': 'X',
             'YRK': 'X',
             'YRM': 'X',
             'YRN': 'X',
             'YRR': 'QR*W',
             'YRS': 'X',
             'YRT': 'HRYC',
             'YRV': 'X',
             'YRW': 'X',
             'YRY': 'HRYC',
             'YSA': 'PRS*',
             'YSB': 'X',
             'YSC': 'PRSC',
             'YSD': 'X',
             'YSG': 'PRSW',
             'YSH': 'X',
             'YSK': 'X',
             'YSM': 'X',
             'YSN': 'X',
             'YSR': 'X',
             'YSS': 'X',
             'YST': 'PRSC',
             'YSV': 'X',
             'YSW': 'X',
             'YSY': 'PRSC',
             'YTA': 'L',
             'YTB': 'LF',
             'YTC': 'LF',
             'YTD': 'LF',
             'YTG': 'L',
             'YTH': 'LF',
             'YTK': 'LF',
             'YTM': 'LF',
             'YTN': 'LF',
             'YTR': 'L',
             'YTS': 'LF',
             'YTT': 'LF',
             'YTV': 'LF',
             'YTW': 'LF',
             'YTY': 'LF',
             'YVA': 'X',
             'YVB': 'X',
             'YVC': 'X',
             'YVD': 'X',
             'YVG': 'X',
             'YVH': 'X',
             'YVK': 'X',
             'YVM': 'X',
             'YVN': 'X',
             'YVR': 'X',
             'YVS': 'X',
             'YVT': 'X',
             'YVV': 'X',
             'YVW': 'X',
             'YVY': 'X',
             'YWA': 'QL*',
             'YWB': 'X',
             'YWC': 'HLYF',
             'YWD': 'X',
             'YWG': 'QL*',
             'YWH': 'X',
             'YWK': 'X',
             'YWM': 'X',
             'YWN': 'X',
             'YWR': 'QL*',
             'YWS': 'X',
             'YWT': 'HLYF',
             'YWV': 'X',
             'YWW': 'X',
             'YWY': 'HLYF',
             'YYA': 'PLS',
             'YYB': 'PLSF',
             'YYC': 'PLSF',
             'YYD': 'PLSF',
             'YYG': 'PLS',
             'YYH': 'PLSF',
             'YYK': 'PLSF',
             'YYM': 'PLSF',
             'YYN': 'PLSF',
             'YYR': 'PLS',
             'YYS': 'PLSF',
             'YYT': 'PLSF',
             'YYV': 'PLSF',
             'YYW': 'PLSF',
             'YYY': 'PLSF'}
        res_triplet_table = self.aligner.generate_table()

        self.assertEqual(exp_triplet_table, res_triplet_table)

    def testEnumerateCodonPossibilities(self):
        # Setting params
        triplet = 'YTD'

        exp_possibilities = ['CTA', 'CTG', 'CTT', 'TTA', 'TTG', 'TTT']
        res_possibilities = self.aligner.enumerate_codon_possibilities(triplet)

        self.assertEqual(exp_possibilities, res_possibilities)

    def testTrimLowQualities(self):
        # Setting params
        codon_list = ['GTC', 'AGT']
        left = 715
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V', 'I'), 37: ('N', 'S', 'N')}
        frameshifts = []
        subtype = ''

        exp_trim_left, exp_right_trim = (0, 0)
        res_trim_left, res_trim_right = \
            self.aligner.trim_low_qualities(codon_list,
                                            left,
                                            first_aa,
                                            last_aa,
                                            mutations,
                                            frameshifts,
                                            gene,
                                            subtype)
        
        self.assertEqual((res_trim_left, res_trim_right),
                         (exp_trim_left, exp_right_trim))

        # Setting params
        codon_list = ['GTC', 'AGT']
        left = 56
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V', 'I'), 37: ('N', 'S', 'N')}
        frameshifts = []
        subtype = 'B'

        exp_trim_left, exp_right_trim = (0, 0)
        res_trim_left, res_trim_right = \
            self.aligner.trim_low_qualities(codon_list,
                                            left,
                                            first_aa,
                                            last_aa,
                                            mutations,
                                            frameshifts,
                                            gene,
                                            subtype)
        
        self.assertEqual((res_trim_left, res_trim_right),
                         (exp_trim_left, exp_right_trim))

        # Setting params
        codon_list = ['GTC', 'AGT']
        left = 155
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V', 'I'), 37: ('N', 'S', 'N')}
        frameshifts = []
        subtype = ''

        exp_trim_left, exp_right_trim = (0, 0)
        res_trim_left, res_trim_right = \
            self.aligner.trim_low_qualities(codon_list,
                                            left,
                                            first_aa,
                                            last_aa,
                                            mutations,
                                            frameshifts,
                                            gene,
                                            subtype)
        
        self.assertEqual((res_trim_left, res_trim_right),
                         (exp_trim_left, exp_right_trim))

    def testIsUnsequenced(self):
        # Setting params
        triplet = 'NNN'
        result = self.aligner.is_unsequenced(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'N-N'
        result = self.aligner.is_unsequenced(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'NNG'
        result = self.aligner.is_unsequenced(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'NTG'
        result = self.aligner.is_unsequenced(triplet)
        self.assertFalse(result)

    def testIsStopCodon(self):
        # Setting params
        triplet = 'TAG'
        result = self.aligner.is_stop_codon(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'TAA'
        result = self.aligner.is_stop_codon(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'TGA'
        result = self.aligner.is_stop_codon(triplet)
        self.assertTrue(result)

        # Setting params
        triplet = 'NNN'
        result = self.aligner.is_stop_codon(triplet)
        self.assertFalse(result)

        # Setting params
        triplet = 'NTG'
        result = self.aligner.is_stop_codon(triplet)
        self.assertFalse(result)

    def testGetHighestMutPrevalence(self):
        # Setting params
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V', 'I'), 37: ('N', 'S', 'N')}
        position = 3
        subtype = ''

        exp_prevalence = 100.0
        res_prevalence = self.aligner.get_highest_mut_prevalance((position, mutations[position]),
                                                                 gene,
                                                                 subtype)

        self.assertEqual(exp_prevalence, res_prevalence)

        # Setting params
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V', 'I'), 37: ('N', 'S', 'N')}
        position = 3
        subtype = ''

        exp_prevalence = 100.0
        res_prevalence = self.aligner.get_highest_mut_prevalance((position, mutations[position]),
                                                                 gene,
                                                                 subtype)

        self.assertEqual(exp_prevalence, res_prevalence)

    def testGetMutPrevalence(self):
        # Setting params
        gene, first_aa, last_aa, first_na, last_na = ('IN', 1, 99, 1, 294)
        mutations = {1: ('F', 'S'), 37: ('N', 'S')}
        position = 1
        cons, aas = mutations[position]
        subtype = 'A'

        expected = 0.0
        result = self.aligner.get_mut_prevalence(position,
                                                 cons,
                                                 aas[0],
                                                 gene,
                                                 subtype)

        self.assertEqual(expected, result)

        # Setting params
        gene, first_aa, last_aa, first_na, last_na = ('PR', 1, 99, 1, 294)
        mutations = {3: ('I', 'V'), 37: ('N', 'S')}
        position = 3
        cons, aas = mutations[position]
        subtype = 'A'

        expected = 0.7
        result = self.aligner.get_mut_prevalence(position,
                                                 cons,
                                                 aas[0],
                                                 gene,
                                                 subtype)

        self.assertEqual(expected, result)

        # Setting params
        gene, first_aa, last_aa, first_na, last_na = ('RT', 1, 99, 1, 294)
        mutations = {4: ('P', 'S'), 37: ('N', 'S')}
        position = 4
        cons, aas = mutations[position]
        subtype = 'F'

        expected = 3.7
        result = self.aligner.get_mut_prevalence(position,
                                                 cons,
                                                 aas[0],
                                                 gene,
                                                 subtype)

        self.assertEqual(expected, result)
