function qf = cui_scheme(qC,qU,qD)

    ac = (qC - qU)/(qD - qU);
    
    if (qD == qU)
        qf = qC;
    else   
        af = [];
        if (0.0 < ac && ac <= 2.0/13.0) 
            af = 3.0*ac;
        elseif (2.0/13.0 < ac && ac <= 4.0/5.0) 
            af = ac*5.0/6.0 + 1/3.0;
        elseif (4.0/5.0 < ac && ac <= 1.0) 
            af = 1.0;
        else
            af = ac;
        end
        qf = af*(qD - qU) + qU;
    end

end

