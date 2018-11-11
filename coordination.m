function coordination( atom,max_coordination_number_plus_one,atomtypes,string_type,Number_of_types )
%coordination analyzer
N=size(atom,2);
bin=zeros(Number_of_types,max_coordination_number_plus_one);
for i=1:N
    for j=1:Number_of_types
        switch atom(i).type
            case atomtypes(j)
                bin(atomtypes(j),atom(i).coordination+1)=bin(atomtypes(j),atom(i).coordination+1)+1;
        end
    end
end
executablestring=['legend({'];
colors=['b' 'r' 'g' 'c' 'm' 'y' 'k'];
figure,
for i=1:Number_of_types
    executablestring=[executablestring char(39) string_type{atomtypes(i)} '=' num2str(round(sum([0:(max_coordination_number_plus_one-1)].*bin(atomtypes(i),:)/sum(bin(atomtypes(i),:))),2)) char(39) ','];
    bar([0:(max_coordination_number_plus_one-1)],bin(atomtypes(i),:)/sum(bin(atomtypes(i),:)),1-i/5,colors(i))
    hold on
end
set(gca,'FontSize',24);
grid on;
hold off;executablestring(end)=[]; 
executablestring = [executablestring '},' char(39) 'FontSize' char(39) ',24);'];
eval(executablestring);
%legend({'O','Si','Zr'},'FontSize',24)
xlabel('Coordination number','fontsize',30)
ylabel(['Probability'],'fontsize',30)
% for i=1:Number_of_types
%     t=text(0,0.9-0.2*i/10,[[string_type(i) '='] num2str(round(sum([0:(max_coordination_number_plus_one-1)].*bin(atomtypes(i),:)/sum(bin(atomtypes(i),:))),2))],'fontsize',24)
%     t.Extent=[0,0.9-0.2*i/10, 0.4, 0.2];
% end
end


