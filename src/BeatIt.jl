module BeatIt

SMILES_CHARS = [' ',
                  '#', '%', '(', ')', '+', '-', '.', '/',
                  '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                  '=', '@',
                  'A', 'B', 'C', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P',
                  'R', 'S', 'T', 'V', 'X', 'Z',
                  '[', '\\', ']',
                  'a', 'b', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'p', 'r', 's',
                  't', 'u']

smi2index = Dict((c,i) for (i,c) in enumerate(SMILES_CHARS))
index2smi = Dict((i,c) for (i,c) in enumerate(SMILES_CHARS))

function smiles_encoder(smiles)
    maxlen=120
    X = zeros(Int, (maxlen, length(SMILES_CHARS)))
    for (i,c) in enumerate(smiles)
        X[i, smi2index[c]] = 1
    end
    return X
end

#function smiles_decoder(X)
#    smi = ""
#    X = maximum(X, dims=1)
#    for i in X
#        smi += index2smi[i]
#    end
#    return smi
#end

end 
