% Resistors/pn/pin: electrical structures
STR{1}='Resistor_ox';               % simple resistor
STR{2}='PN_ox';                     % PN junction, with oxide aperture
STR{3}='PIN_ox';                    % PIN junction, with oxide aperture
STR{4}='Resistor_noOX_simple';      % simple resistor

% JSTQE 
STR{6}='MarkusN_FINALEdisFitto';    % Markus, for JSTQE article!
STR{7}='MarkusN_TopLAqw2';          % 1D version, some adjustments (used in D1ANA at the beginning)
STR{8}='MarkusN_FINALE_LgDBR';      % Markus, for JSTQE article; Lg in DBR!
STR{9}='MarkusN_FINALE_LgDBR_fixed'; % Markus, for JSTQE article; Lg in DBR!

% 1TJ etched (lithographic, NO OX!) - double n-DBR
% "Litho 1"
STR{10}='MarkusN_BTJetch_opt_VERDE';    % lithographic. TJ (only in VELM, no mireq)
STR{11}='MarkusN_BTJetch_VERDE_LgDBR';  % lithographic. TJ (only in VELM, Lg DBR)
STR{12}='MarkusN_BTJetch_VERDE';
STR{13}='MarkusN_BTJetch_opt_mireq';    % lithographic. TJ (only in VELM, w/out top mireq)

% "Litho 2"
STR{15}='MarkusN_BTJetch_DD';           % lithographic. TJ (both VELM and DD)
STR{16}='MarkusN_BTJetch_DD_ROSSO_LgDBR';     % lithographic. TJ (both VELM and DD, Lg DBR)
STR{17}='MarkusN_BTJetch_DD_ROSSO';     % lithographic. TJ ("25", both VELM and DD; fewer mesh points)
STR{18}='MarkusN_BTJetch_DD_ROSSO_LgDBR_singleMode';     % lithographic. TJ (both VELM and DD, Lg DBR)
STR{19}='MarkusN_BTJetched25eq';

% 1TJ (infinite + OX in different positions - (the oxide in radially graded as in JSTQE Markus!)) 
STR{20}='MarkusN_TJ_oxAbove';
STR{21}='MarkusN_TJ_oxBelow';
STR{22}='MarkusN_TJ_oxBoth';

% 1TJ (Lg in DBR!, as in MarkusN_FINALE_LgDBR.str)
STR{25}='MarkusN_TJ_oxAbove_LgDBR';
STR{26}='MarkusN_TJ_oxBelow_LgDBR';
STR{27}='MarkusN_TJ_oxAbove_LgDBR_fixed';
STR{28}='MarkusN_TJ_oxBelow_LgDBR_fixed';

STR{29}='MarkusN_TJ_oxAbove_LgDBR_oxAnti';

% 1TJ (with p-DBR and OX aperture: only current recycling)
STR{30}='MarkusN_TJ_pDBR_2AR';      % maybe from MarkusN_2TJ_inf_BELOW.str
STR{31}='MarkusN_TJ_pDBR_2AR_alternative';      % maybe from MarkusN_2TJ_inf_BELOW.str

% 1TJ (radially limited by OX, like MarkusN_BTJetch_opt (w/out central step in VELM) - SIMPLIFICATION!)
STR{32}='MarkusN_TJ_radialOx';      % obtained FROM 'MarkusN_BTJetch_opt.str'

% 2TJ
STR{35}='MarkusN_2TJ';                      % oxide beside each TJ

STR{36}='MarkusN_2TJ_inf_OXbelowAR';        % radially unlimited TJs, oxides BELOW each TJ+AR
STR{37}='MarkusN_2TJ_inf+LgDBR_OXbelowAR';  % radially unlimited TJs, oxides BELOW each TJ+AR; Lg mirrors!

STR{38}='MarkusN_2TJ_inf_OXaboveAR';        % radially unlimited TJs, oxides ABOVE each TJ+AR
STR{39}='MarkusN_2TJ_inf+LgDBR_OXaboveAR';  % radially unlimited TJs, oxides BELOW each TJ+AR; Lg mirrors!

% 3TJ
STR{40}='MarkusN_3TJ';                      % oxide beside every TJ 
STR{41}='MarkusN_3TJ_inf';                  % radially unlimited TJs, 3 oxides BELOW each TJ
STR{42}='MarkusN_3TJ_inf+LgDBR_OXaboveAR';  % radially unlimited TJs, oxides ABOVE each TJ+AR

% Stephan structures (GR and SR)
STR{50}='Stephan_GRok';     % Metal
STR{51}='Stephan_GR1';     % Torrelli GR

% Julian TJ-VCSEL
STR{60}='Julian_TJ_2AR_pDBR';
STR{61}='Julian_TJ_2AR_pDBR_fake';


