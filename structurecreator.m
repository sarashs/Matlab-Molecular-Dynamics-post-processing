%clear all
%% Global variables 
% Number_of_atoms, Number_of_types
%% End of global
while 1
    UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');
    set(0, 'DefaultUIControlFontSize', 18);
    choice = menu('Choose a calculation','Load data','RDF','Neighbor list calculation','Angle analysis','Coordination','Save workspace','Percentage of atoms','Charge distribution/Dipole analysis','Clear workspace','Exit');
    
    switch choice
        case 1
            set(0, 'DefaultUIControlFontSize', 18);
            choice2 = menu('Data input options','Load workspace from previous calculation','Load your new extended xyz data file');
            switch choice2
                case 1
                    [filename, pathname] = uigetfile('*.mat', 'Select a MATLAB .mat file.');
                    if isequal(filename,0)
                        disp('User selected Cancel')
                    else
                        disp(['User selected ', fullfile(pathname, filename)])
                    end
                    load(fullfile(pathname, filename))
                case 2
                    %% Import from Ovito output
                    [ importedxyz, filename, pathname ]  = importovito;
                    %%
                    Number_of_atoms=importedxyz(1,1);
                    maxx=importedxyz(2,2); maxy=importedxyz(2,6); maxz=importedxyz(2,10);
                    for i=3:(Number_of_atoms+2)
                        
                        atom(i-2).id=importedxyz(i,1);
                        atom(i-2).type=importedxyz(i,2);
                        atom(i-2).xyz=importedxyz(i,4:6);
                        atom(i-2).charge=importedxyz(i,3);
                        %%%%%%%%premake a structure
                    end
                    clearvars i A importedxyz
                    savename=[filename(1:end-4) '.mat'];
                    %save(fullfile(pathname, savename),'atom')
                    
                    %%%%NOTE the given vector is from the atom to its' neighbors
                    %clearvars namenumber
                    %         for i=1:Number_of_types
                    %             prompt = ['Numerical value for type ' int2str(i) 'st is:'];
                    %             atomtypes(i)=input(prompt)
                    %         end
                    %% Input the atom types
                    T=struct2table(atom);
                    atomtypes=unique(T.type);
                    clearvars T;
                    string_type=[]; Number_of_types=size(atomtypes,1);
                    executablestring = ['prompt = {'];
                    for i=1: Number_of_types
                        executablestring=[executablestring char(39) 'Type ' num2str(atomtypes(i)) ' element:' char(39) ','];
                    end
                    executablestring(end)=[];
                    executablestring=[executablestring '};'];
                    eval(executablestring);
                    dlg_title = 'Element/type characterization';
                    line_size=repmat([1 20],Number_of_types);line_size=line_size(:,1:2);
                    answer = inputdlg(prompt,dlg_title,line_size);
                    string_type={};
                    for i=1:1:Number_of_types
                        string_type{atomtypes(i)}=answer{i};
                    end
                    clearvars answer choice dlg_title i executablestring line_size prompt UIControl_FontSize_bak
            end
            %% RDF Calculation
        case 2
            %% RDF max radius etc input
            dlg_title = 'RDF Essentials'; prompt={'Maximum radius (1/2 of the smallest box dimension):','Step (0.2 recommended)'};
            line_size=repmat([1 140],2);line_size=line_size(:,1:2);
            set(0, 'DefaultUIControlFontSize', 18);
            answer = inputdlg(prompt,dlg_title,line_size);
            Maximum_RDF_Radius= str2double(answer{1}); Step= str2double(answer{2});
            clearvars prompt dlg_title line_size
            %% which RDF input
            String_of_calculations=[];
            executablestring = ['List_String_of_calculations = {'];
            for i=1:Number_of_types
                for j=i:Number_of_types
                    executablestring = [executablestring char(39) string_type{atomtypes(i)} '-' string_type{atomtypes(j)} char(39) ','];
                    String_of_calculations(end+1,:)=[atomtypes(i),atomtypes(j)];
                end
            end
            executablestring(end) = []; executablestring = [executablestring '};'];
            eval(executablestring);
            [set_of_calculations,v] = listdlg('PromptString','Select the RDFs to be calculated:',...
                'SelectionMode','multiple',...
                'ListString',List_String_of_calculations);
            clearvars  i j List_String_of_calculations executablestring
            %% RDF computations
            radius=0:Step:Maximum_RDF_Radius;Number_of_atoms=size(atom,2);S=size(radius,2);
            Volume=maxx*maxy*maxz;
            for k=1:size(set_of_calculations,2)
                particle_frequency=0;
                eval(['RDF' string_type{String_of_calculations(set_of_calculations(k),1)} string_type{String_of_calculations(set_of_calculations(k),2)} '=zeros(S,1);']);
                for i=1:Number_of_atoms
                    if atom(i).type==String_of_calculations(set_of_calculations(k),1)
                        particle_frequency=particle_frequency+1;
                        for  j=1:Number_of_atoms
                            if atom(j).type==String_of_calculations(set_of_calculations(k),2)
                                [ x, y, z ] = subtract( atom(j).xyz, atom(i).xyz, maxx, maxy, maxz ); A=sqrt([ x, y, z ]*[ x, y, z ]');
                                for r=1:S
                                    if A > radius(r) && A <= (radius(r)+Step) && i~=j
                                        eval(['RDF' string_type{String_of_calculations(set_of_calculations(k),1)} string_type{String_of_calculations(set_of_calculations(k),2)} '(r)=RDF' string_type{String_of_calculations(set_of_calculations(k),1)} string_type{String_of_calculations(set_of_calculations(k),2)} '(r)+1/(Step*4*pi*radius(r)^2);']);
                                    end
                                end
                            end
                        end
                    end
                end
                eval(['RDF' string_type{String_of_calculations(set_of_calculations(k),1)} string_type{String_of_calculations(set_of_calculations(k),2)} '=RDF' string_type{String_of_calculations(set_of_calculations(k),1)} string_type{String_of_calculations(set_of_calculations(k),2)} '*particle_frequency/Volume;']);
            end
            clearvars i j k Step r N A x y z
            %% RDF Display
            for p=1:size(set_of_calculations,2)
                eval(['Display_var=RDF' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} ';']);
                figure, plot(radius',Display_var,'-.','linewidth',3);
                xlim([0 Maximum_RDF_Radius]);
                %xticks([0 30 60 90 120 150 180]);
                grid on;
                set(gca,'FontSize',24);
                xlabel('r','fontsize',30);
                ylabel('g_{ab}(r)','fontsize',30);
                legend({[string_type{String_of_calculations(set_of_calculations(p),1)} '-' string_type{String_of_calculations(set_of_calculations(p),2)} ]},'FontSize',24);
            end
            %clearvars Maximum_RDF_Radius answer p Display_var Step line_size prompt dlg_title set_of_calculations v String_of_calculations
            %% Neighbor list calculation
        case 3
            executablestring = ['prompt = {']; num_lines=0;
            executablestring2 = ['defaultans = {'];
            bond_distances=zeros(Number_of_types*(Number_of_types-1),Number_of_types*(Number_of_types-1));
            bonds=zeros(Number_of_types*(Number_of_types-1),3);
            for i=1: Number_of_types
                for j=1: Number_of_types
                    if j>=i
                        num_lines=num_lines+1;
                        bonds(num_lines,1:2)=[atomtypes(i) atomtypes(j)];
                        executablestring=[executablestring char(39) string_type{atomtypes(i)} '-' string_type{atomtypes(j)} ' maximum bond distance:' char(39) ','];
                        executablestring2 = [ executablestring2 char(39) '0' char(39) ','];
                    end
                end
            end
            executablestring(end)=[]; executablestring=[executablestring '};'];
            executablestring2(end) = []; executablestring2 = [executablestring2 '};'];
            eval(executablestring);
            eval(executablestring2);
            dlg_title = 'Bond Distance Input';
            line_size=repmat([1 40],num_lines);line_size=line_size(:,1:2);
            answer = inputdlg(prompt,dlg_title,line_size,defaultans);
            bonds(:,3)= str2double(answer);
            for i=1:num_lines
                bond_distances(bonds(i,1),bonds(i,2))=bonds(i,3);
            end
            bond_temp=[];
            for i=1:num_lines
                if bonds(i,3)~=0
                    bond_temp(end+1,:)=bonds(i,:);
                end
            end
            clearvars bonds
            bonds=bond_temp;
            bonds_size=size(bonds,1);
            clearvars bond_temp answer prompt line_size dlg_title executablestring
            neighbor(1)=struct('type',0,'id',0,'vector',zeros(3,1));scount=24; % we assume no atom has more than 24 neighbors
            for neigh_count=2:scount
                neighbor(neigh_count)=neighbor(1);
            end
            k=0; Number_of_atoms=size(atom,2);
            tic
            for i=1:Number_of_atoms
                atom(i).neighbors=neighbor; atom(i).coordination=0;
            end
            h = waitbar(0,'Initialization ...');
            for i=1:Number_of_atoms
                for j=1:Number_of_atoms
                    for k=1:bonds_size
                        if atom(i).type == bonds(k,1) && i~=j && atom(j).type==bonds(k,2)
                            [ x, y, z ] = subtract( atom(j).xyz, atom(i).xyz, maxx, maxy, maxz ); A=sqrt([ x, y, z ]*[ x, y, z ]');
                            if A < bond_distances(bonds(k,1),bonds(k,2))
                                atom(i).coordination=atom(i).coordination+1;
                                atom(i).neighbors(atom(i).coordination).id=atom(j).id;
                                atom(i).neighbors(atom(i).coordination).type=atom(j).type;
                                atom(i).neighbors(atom(i).coordination).vector=[x, y, z ];
                                
                                if bonds(k,1)~=bonds(k,2)
                                    atom(j).coordination=atom(j).coordination+1;
                                    atom(j).neighbors(atom(j).coordination).id=atom(i).id;
                                    atom(j).neighbors(atom(j).coordination).type=atom(i).type;
                                    atom(j).neighbors(atom(j).coordination).vector=-[x, y, z ];
                                    waitbar(i/Number_of_atoms,h,sprintf('%f%% along...',i/Number_of_atoms*100))
                                end
                            end
                        end
                    end
                end
            end
            delete(h)
            toc
            clearvars neighbor scount neigh_count bonds_size
            %save(fullfile(pathname, savename),'atom');
            %% Angle
        case 4
            %dlg box to decide upon the calculation
            String_of_calculations=[];
            executablestring = ['List_String_of_calculations = {'];
            for i=1:Number_of_types
                for j=1:Number_of_types
                    for k=1:Number_of_types
                        executablestring = [executablestring char(39) string_type{atomtypes(i)} '-' string_type{atomtypes(j)} '-' string_type{atomtypes(k)} char(39) ','];
                        String_of_calculations(end+1,:)=[atomtypes(i),atomtypes(j),atomtypes(k)];
                    end
                end
            end
            executablestring(end) = []; executablestring = [executablestring '};'];
            eval(executablestring);
            [set_of_calculations,v] = listdlg('PromptString','Select the angles to be calculated:',...
                'SelectionMode','multiple',...
                'ListString',List_String_of_calculations);
            %set_of_calculations is the set of calculations that should be performed from String_of_calculations
            %%
            T = struct2table(atom);
            [occurances,types]=hist(T.type,unique(T.type));
            S=size(atom,2); max_coordination_number=max(T.coordination);
            Number_of_atoms = max_coordination_number*S;
            Number_of_angle_calculations=size(set_of_calculations,2);
            for i=1:Number_of_angle_calculations
                eval(['data' string_type{String_of_calculations(set_of_calculations(i),1)} string_type{String_of_calculations(set_of_calculations(i),2)} string_type{String_of_calculations(set_of_calculations(i),3)} '=zeros(Number_of_atoms,1);']);
                eval(['pntr' string_type{String_of_calculations(set_of_calculations(i),1)} string_type{String_of_calculations(set_of_calculations(i),2)} string_type{String_of_calculations(set_of_calculations(i),3)} '=1;']);
            end
            for i=1:Number_of_types
                eval(['IDcheck' num2str(atomtypes(i)) '=zeros(S,1);']);
            end
            for i=1:S
                atom(i).angle=zeros(max_coordination_number,1);
            end
            for p=1:Number_of_angle_calculations
                for i=1:S
                    if atom(i).coordination~=0 && atom(i).type==String_of_calculations(set_of_calculations(p),1) && eval(['IDcheck' num2str(String_of_calculations(set_of_calculations(p),1)) '(atom(i).id)~=1'])
                        for j=1:atom(i).coordination
                            if atom(atom(i).neighbors(j).id).coordination>1 && atom(i).neighbors(j).type==String_of_calculations(set_of_calculations(p),2)
                                for k=1:atom(atom(i).neighbors(j).id).coordination
                                    if atom(atom(i).neighbors(j).id).neighbors(k).type==String_of_calculations(set_of_calculations(p),3) && atom(atom(i).neighbors(j).id).neighbors(k).id ~= atom(i).id && eval(['IDcheck' num2str(String_of_calculations(set_of_calculations(p),3)) '((atom(atom(i).neighbors(j).id).neighbors(k).id))~=1'])
                                        ANGL  = angles(-atom(i).neighbors(j).vector,atom(atom(i).neighbors(j).id).neighbors(k).vector);
                                        eval(['pntr=pntr' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} string_type{String_of_calculations(set_of_calculations(p),3)} ';']);
                                        eval(['data' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} string_type{String_of_calculations(set_of_calculations(p),3)} '(pntr)=ANGL;']);
                                        eval(['pntr' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} string_type{String_of_calculations(set_of_calculations(p),3)} '=pntr+1;']);
                                        eval(['IDcheck' num2str(String_of_calculations(set_of_calculations(p),3)) '((atom(atom(i).neighbors(j).id).neighbors(k).id))=1;']);
                                        eval(['IDcheck' num2str(String_of_calculations(set_of_calculations(p),1)) '(atom(i).id)=1;']);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            for p=1:Number_of_angle_calculations
                eval(['angledata=data' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} string_type{String_of_calculations(set_of_calculations(p),3)} ';']);
                angledata=angledata(1:find(angledata, 1, 'last' ));
                figure, hist(angledata,15);
                xlim([0 180]);
                xticks([0 30 60 90 120 150 180]);
                grid on;
                set(gca,'FontSize',24);
                eval(['data' string_type{String_of_calculations(set_of_calculations(p),1)} string_type{String_of_calculations(set_of_calculations(p),2)} string_type{String_of_calculations(set_of_calculations(p),3)} '=angledata;']);
                xlabel('Angle','fontsize',30);
                ylabel('Frequency','fontsize',30);
                legend({[string_type{String_of_calculations(set_of_calculations(p),1)} '-' string_type{String_of_calculations(set_of_calculations(p),2)} '-' string_type{String_of_calculations(set_of_calculations(p),3)}]},'FontSize',24)
            end
            clearvars occurances types T S angledata i j k List_String_of_calculations String_of_calculations v
            %% Coordination
        case 5
            coordination(atom,9,atomtypes,string_type,Number_of_types);
            %% Save data
        case 6
            save(fullfile(pathname, savename));
            %% Percentage of types
        case 7
            T = struct2table(atom);
            [occurances,types]=hist(T.type,unique(T.type));
            occurances=occurances/sum(occurances);
            fileID = fopen('percentage.txt','a'); 
            executablestring = [char(39)];
            for i=1:size(occurances,2)
                executablestring=[executablestring '%5.3f '];
            end
            executablestring =[executablestring '\n' char(39) ',types,occurances);'];
            fprintf(fileID,[filename '\n']);
            filename
            eval(['fprintf(fileID,' executablestring]);
            occurances
            %fprintf(fileID,'%4.2f %4.2f %4.2f\n',types,occurances);
            fclose(fileID);
            clearvars T %executablestring
            %% Charge distribution
        case 8
            direction = menu('Choose a Direction','X','Y','Z');
            prompt = 'Window size in angstrom? ';
            bin_size = input(prompt);
            [charge, Voltage, direction_axis, odensity]=chargedistribution(atom, direction,bin_size,maxx,maxy,maxz);
        case 9
            clear all
            %    for i=1:Number_of_types
            %        eval(['clearvars IDcheck' num2str(atomtypes(i)) ';']);
            %    end
            %clearvars N filename choice atomtypes maxx maxy maxz savename delimiter UIControl_FontSize_bak Volume set_of_calculations pntr pathname particle_frequency p Number_of_types i l Display_var executablestring
        case 10
            break
    end
end