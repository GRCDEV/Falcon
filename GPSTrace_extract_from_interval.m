function trace2 = GPSTrace_extract_from_interval(trace,t1,t2)
    % Extract GPS trace from t1 to t2
    for i = 1:length(trace)
        if trace(i,2) >= t1
            break;
        end
    end

    iStart = i;

    for i = iStart:length(trace)
        if trace(i,2) >= t2
            break;
        end
    end
    
    iEnd = i;

    trace2 = trace(iStart:iEnd,:);
end
