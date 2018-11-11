function [charge_distribution, Voltage, direction_axis, odensity]=chargedistribution(atom,direction,bin_size,maxx,maxy,maxz);
switch direction
    case 1
        direction_axis=0:(maxx/round(maxx/bin_size)):maxx;
        area=maxy*maxz;
    case 2
        direction_axis=0:(maxy/round(maxy/bin_size)):maxy;
        area=maxx*maxz;
    case 3
        direction_axis=0:(maxz/round(maxz/bin_size)):maxz;
        area=maxy*maxx;
end

ii=size(direction_axis,2); mid=(direction_axis(2)-direction_axis(1)/2);
charge_distribution=zeros(ii-1,1);
odensity=zeros(ii-1,1);
%charge_distribution_freq=zeros(ii-1,1);
Number_of_atoms=size(atom,2);
Voltage=0;
for j=1:Number_of_atoms
        if atom(j).xyz(direction)>direction_axis(end)
            temp=atom(j).xyz(direction)-direction_axis(end);
        elseif atom(j).xyz(direction)<0
            temp=atom(j).xyz(direction)+direction_axis(end);
        else
            temp=atom(j).xyz(direction);
        end
    for i=2:ii
        if temp<(direction_axis(i)) && temp>=(direction_axis(i-1))
            charge_distribution(i-1)=charge_distribution(i-1)+atom(j).charge;
            if atom(j).type==1
                odensity(i-1)=odensity(i-1)+1;
            end
 %           charge_distribution_freq(i-1)=charge_distribution_freq(i-1)+1;
            break
        end
    end
end
%charge_distribution=charge_distribution./charge_distribution_freq;
%Voltage=;
direction_axis=direction_axis+mid;
direction_axis(end)=[];
%odensity=1E24*odensity/(mid*2*area); %atoms per cm^3
charge_distribution=1E24*1.602E-19*charge_distribution/(mid*2*area); %C per cm^3
%figure,plot(direction_axis, odensity);
figure, plot(direction_axis,charge_distribution);
%figure, plot(direction_axis,charge_distribution./charge_distribution_freq);
end

