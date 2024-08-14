function w = nuttallwin(N)
persistent nutTable

if isempty(nutTable),
    load nuttwintable
end

w = nutTable{N};
