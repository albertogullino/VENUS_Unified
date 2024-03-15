clear
close all
clc

load Data_VENUSThermalSources-2mA.mat

NonMatrixData = [I, V, Pst, PTherm, PJoule, Pfca, Pnr, Prad, Pcap, DeltaTmax, DeltaTmax_Joule, DeltaTmax_FCA, DeltaTmax_nr, DeltaTmax_rad, DeltaTmax_Ccap];     

xlswrite('Data_VENUSThermalSources-Curr.xlsx',NonMatrixData)