using PyCall
using Statistics
using PyPlot
np = pyimport("numpy")

#Disordered initial condition
function initialization(Ly, Lx, rho)
        #config = zeros(Ly, Lx)
        config = fill(-1, (Ly, Lx))
        pcount = 0
        numOfPart = convert(Int64, floor(Lx*Ly*rho))
        while (pcount < numOfPart)
            i = rand(1:Ly)
            j = rand(1:Lx)
            if (config[i, j] == -1)
                pcount += 1
                config[i, j] = 1
            end
        end
        checkNum(config, Ly, Lx)
        return config
end

#Horizontally ordered initial config
function initialization2(Ly, Lx, rho)
        #config = zeros(Ly, Lx)
        config = fill(-1, (Ly, Lx))
        numOfOneType = convert(Int64, floor(Lx*Ly*rho/2))
        counter = numOfOneType*2
        i, j = 1, 1
        while counter > 0
                if j == Lx
                    config[i, j] = 1
                    i += 1
                    j = 1
                    counter -= 1
                else
                    config[i, j] = 1
                    j += 1
                    counter -= 1
                end
        end
        return config
end

#Vertically ordered initial config with gap
function initialization3(Ly, Lx, rho)
        #config = zeros(Ly, Lx)
        config = fill(-1, (Ly, Lx))
        numOfOneType = convert(Int64, floor(Lx*Ly*rho/2))
        counter = numOfOneType*2
        i, j = 1, 1
        while counter > 0
            if i == Lx
                config[i, j] = 1
                j += 1
                i = 1
                counter -= 1
            else
                config[i, j] = 1
                i += 1
                counter -= 1
            end
        end
        return config
end

function moveCheck(config, Ly, Lx, Efield, beta)
    for i = 1:Ly*Lx
        a = rand(1: Ly)
        b = rand(1: Lx)
        if config[a, b] == 1
            dirc = rand(1: 4)
            if dirc == 1
                if config[a, 1+mod(b,Lx)] == -1
                    config =  mcmoveph(config, Ly, Lx, a, b, 1, Efield, beta)
                end
            end
            if dirc == 2
                if config[1+mod(a,Ly),b] == -1
                    config =  mcmoveph(config, Ly, Lx, a, b, 2, Efield, beta)
                end
            end
            if dirc == 3
                if config[a,1+mod(b-2,Lx)] == -1
                    config =  mcmoveph(config, Ly, Lx, a, b, 3, Efield, beta)
                end
            end
            if dirc == 4
                if config[1+mod(a-2,Ly),b] == -1
                    config =  mcmoveph(config, Ly, Lx, a, b, 4, Efield, beta)
                end
            end
        end
    end
    return config

end

function mcmoveph(config, Ly, Lx, a, b, dirc, Efield, beta)
    E = Efield
    s =  config[a, b]
    above = config[1+mod(a,Ly),b]
    right = config[a, 1+mod(b,Lx)]
    below = config[1+mod(a-2,Ly),b]
    left = config[a,1+mod(b-2,Lx)]
    """here are the definitions of moves"""
    if s == 1
        if dirc == 1
            nb1 = config[a,1+mod(b-2,Lx)]+config[1+mod(a,Ly),b]+config[1+mod(a-2,Ly),b] 
            nb2 = config[1+mod(a,Ly),1+mod(b,Lx)]+config[1+mod(a-2,Ly),1+mod(b,Lx)]+config[a,1+mod(b+1,Lx)]
            cost = 2*(s*nb1 + right*nb2)      
            if rand()< min(1, exp(-beta*cost))
                s, right = right, s
            end
        end
        if dirc == 2
            nb1 = config[a,1+mod(b-2,Lx)]+config[a, 1+mod(b,Lx)]+config[1+mod(a-2,Ly),b]
            nb2 = config[1+mod(a+1,Ly),b]+config[1+mod(a,Ly),1+mod(b,Lx)]+config[1+mod(a,Ly),1+mod(b-2,Lx)]
            cost = 2*(s*nb1 + above*nb2)
            if rand()< min(1, exp(-beta*(cost-E)))
                s, above = above, s
            end
        end
        if dirc == 3
            nb1 = config[1+mod(a,Ly),b]+config[a, 1+mod(b,Lx)]+config[1+mod(a-2,Ly),b]
            nb2 = config[1+mod(a,Ly),1+mod(b-2,Lx)]+config[a,1+mod(b-3,Lx)]+config[1+mod(a-2,Ly),1+mod(b-2,Lx)]
            cost = 2*(s*nb1 + left*nb2)
            if rand()< min(1, exp(-beta*cost))
                s, left = left, s
            end
        end
        if dirc == 4
            nb1 = config[1+mod(a,Ly),b]+config[a,1+mod(b-2,Lx)]+config[a, 1+mod(b,Lx)]
            nb2 = config[1+mod(a-3,Ly),b]+config[1+mod(a-2,Ly),1+mod(b,Lx)]+config[1+mod(a-2,Ly),1+mod(b-2,Lx)]
            cost = 2*(s*nb1 + below*nb2)
            if rand()< min(1, exp(-beta*(cost+E)))
                s, below = below, s
            end
        end
    end
    config[a, b] = s
    config[1+mod(a,Ly),b] = above 
    config[a, 1+mod(b,Lx)] = right
    config[1+mod(a-2,Ly),b] = below
    config[a,1+mod(b-2,Lx)] = left
    return config
end


function orderParaPsiH(config, Ly, Lx)
        Psi = 0
        for j = 1:Ly
            for i = 1:Lx
                if config[i,j] == 1
                    Psi += exp(2*pi*1im*j/Ly)
                elseif config[i,j] == 0
                    Psi += exp(2*pi*1im*j/Ly)*(-1)
                end
            end
        end
        Psi = abs(Psi)
        Psi = Psi/(Lx*Ly)
        return Psi
end

function orderParaPsiV(config, Ly, Lx)
        Psi = 0
        for j = 1:Ly
            for i = 1:Lx
                if config[i,j] == 1
                    Psi += exp(2*pi*1im*i/Ly)
                elseif config[i,j] == 0
                    Psi += exp(2*pi*1im*i/Ly)*(-1)
                end
            end
        end
        Psi = abs(Psi)
        Psi = Psi/(Lx*Ly)
        return Psi
end

function checkNum(config, Ly, Lx)
    pcount = 0
    hcount = 0
    for i = 1:Ly
        for j = 1:Lx
            if config[i, j] == 1
                pcount += 1
            end
            if config[i, j] == -1
                hcount += 1
            end
        end
    end
end

function configPlot(f, config, Ly, Lx, n_, i)
    sp =  f.add_subplot(1, 6, n_)    
    pcolor(config, cmap=plt.cm.Blues, edgecolors="black") 
    StepString = string("t = ", i)
    plt.title(StepString, fontsize=26)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
end

function simulate(Ly, Lx, rho, Efield, temp, steps)
    config = initialization(Ly, Lx, rho)
    beta = 1.0/temp
    listOfPsiH, listOfPsiV = zeros(0), zeros(0)
    
    f = figure(figsize=(65, 10), dpi=80);    
    configPlot(f, config, Ly, Lx, 1, 0);
    
    MCTimeStep = steps
    for i = 1:MCTimeStep
        config = moveCheck(config, Ly, Lx, Efield, beta) 
        append!(listOfPsiH, orderParaPsiH(config, Ly, Lx))
        append!(listOfPsiV, orderParaPsiV(config, Ly, Lx))
        
        if i == Int(steps*0.0001)       
            configPlot(f, config, Ly, Lx, 2, i)
        end        
        if i == Int(steps*0.001)       
            configPlot(f, config, Ly, Lx, 3, i)
        end      
        if i == Int(steps*0.01)       
            configPlot(f, config, Ly, Lx, 4, i)
        end  
        if i == Int(steps*0.1)       
            configPlot(f, config, Ly, Lx, 5, i)
        end  
        if i == Int(steps)       
            configPlot(f, config, Ly, Lx, 6, i)
        end  
        
    end
    #FigPara = string("E",Efield,"Rho",rho,"temp",temp,"Size",Ly,Lx,".pdf")
    #plt.savefig(FigPara)
    return listOfPsiH, listOfPsiV
    #return config
end


function OvsTemp(Ly, Lx, rho, Efield, minTemp, maxTemp, steps)
    listOfOrderPsiH = zeros(0)
    listOfVarPsiH = zeros(0)
    OPPsiH = zeros(0)
    listOfOrderPsiV = zeros(0)
    listOfVarPsiV = zeros(0)
    OPPsiV = zeros(0)
    
    Temp = range(minTemp; stop=maxTemp, length=10)
    for temp in Temp
        println("Now the temperature value is given by ", temp)
        OPPsiH, OPPsiV = simulate(Ly, Lx, rho, Efield, temp, steps)
        append!(listOfOrderPsiH, Statistics.mean(OPPsiH[500000:end]))
        append!(listOfVarPsiH, Statistics.var(OPPsiH[500000:end]))
        append!(listOfOrderPsiV, Statistics.mean(OPPsiV[500000:end]))
        append!(listOfVarPsiV, Statistics.var(OPPsiV[500000:end]))
    end

    psiH = np.array(listOfOrderPsiH)
    psiH = np.transpose(psiH)
        
    psiHVar = np.array(listOfVarPsiH)
    psiHVar = np.transpose(psiHVar)

        
    
    psiV = np.array(listOfOrderPsiV)
    psiV = np.transpose(psiV)
        
    psiVVar = np.array(listOfVarPsiV)
    psiVVar = np.transpose(psiVVar)
        
    
    return psiH, psiHVar, psiV, psiVVar
end

psiH, psiHVar, psiV, psiVVar = OvsTemp(30, 30, 0.5, 5, 1.8, 5, 1000000)

minTemp = 1.5
maxTemp = 4.0
Temp = range(minTemp; stop=maxTemp, length=10)

plt.plot(Temp, psiV)
plt.plot(Temp, psiH)

plt.plot(Temp, psiHVar)
plt.plot(Temp, psiVVar)