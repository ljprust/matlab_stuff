clear,clc;
results=zeros(1,3);
for i=1:100000
    attacker=randi(6,1,3);
    defender=randi(6,1,2);
    attackermax=max(attacker);
    defendermax=max(defender);
    defendermin=min(defender);
    for j=1:3
        if attacker(j)==attackermax
            attacker(j)=0;
            break
        end
    end
    attacker2nd=max(attacker);
    if attackermax>defendermax && attacker2nd>defendermin
        results(3)=results(3)+1;
    elseif defendermax>=attackermax && defendermin>=attacker2nd
        results(1)=results(1)+1;
    else
        results(2)=results(2)+1;
    end
end
results=results/sum(results)*100
bar(results);