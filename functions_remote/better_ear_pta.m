% order of better_ear corresponds to order of id_patient_list 
% T_pta: table where ptaL and ptaR are stored (and id_patient) - doesnt not
% necessarily be the same as id_patient_list -> list could be a subset for
% example 
% 
% MB 30/05/23
function [better_ear] = better_ear_pta(T_pta,id_patient_list)

% choose better ear for all patients -> to filter data (conditions) for
% individual patients (1: right, 2: left)
better_ear = nan(length(id_patient_list),1); 
for n = 1:length(id_patient_list)
    if T_pta.ptaL(T_pta.id_patient == id_patient_list(n)) >= T_pta.ptaR(T_pta.id_patient == id_patient_list(n))
        better_ear(n,1) = 1; 
    elseif T_pta.ptaL(T_pta.id_patient == id_patient_list(n)) < T_pta.ptaR(T_pta.id_patient == id_patient_list(n))
        better_ear(n,1) = 2; 
    elseif isnan(T_pta.ptaL(T_pta.id_patient == id_patient_list(n))) && ~isnan(T_pta.ptaR(T_pta.id_patient == id_patient_list(n)))
        better_ear(n,1) = 1; 
    elseif ~isnan(T_pta.ptaL(T_pta.id_patient == id_patient_list(n))) && isnan(T_pta.ptaR(T_pta.id_patient == id_patient_list(n)))
        better_ear(n,1) = 2; 
    end
end