function [Data, Subject, Pos, str_ind, pr_ind, tb, Subject_name] = capture_data(GEO_number, Data_GEO,gid,titles,Info,PInfo,geoStruct)

  %Ask which condition to analyze
  tb = tabulate(titles);
  dis = strcat(cellstr(arrayfun(@num2str, 1:length(tb(:,1)), 'UniformOutput', false))',' : ',tb(:,1));
  display(dis);    
	
  display(sprintf('\n\nAll the samples associated to %s are listed above. You must now specify the samples that you want to include in the analysis.\n', GEO_number));
  display(sprintf('You must do this by entering a term or string that is common ONLY to the desired samples from the list above.\n'));
  prompt = 'Which samples would you like to include in the analysis? (enter the common string) ';
  str_ind = input(prompt);
  pr_ind = ExtractSubjects(str_ind, dis);
  display(sprintf('\nThe samples that match your search terms are listed below.\n'));
  display(dis(pr_ind));

  %Get Data for that condtion
  Data     = Data_GEO(:,pr_ind);  
  
  display(sprintf('\nThe list below shows the information available for each one of the selected samples.\n\n'));
    
  %Display the Characteristics of the GEO series
  display(strcat(cellstr(arrayfun(@num2str, 1:length({Info{:,pr_ind(1)}}), 'UniformOutput', false))',' : ',{Info{:,pr_ind(1)}}'));

  %Find out where the time is
  prompt = 'Which item in the list corresponds to the time point? (Enter item number or -1 if none): ';
  tm_ind = input(prompt);
  
  if tm_ind > 0
    % Traditional case: time points can be read from characteristics matrix.
    format bank;
    Pos    = {Info{tm_ind,pr_ind}};
    if ~isempty(cell2mat(strfind(Pos,'Baseline')))
      Pos = strrep(Pos, 'Baseline', '0');
    end
    Pos    = cell2mat(cellfun(@(x) str2num(char(regexp(x,'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match'))), Pos, 'UniformOutput', false));

    sane_check = 0;
    while sane_check == 0
      display(Pos');
      prompt = '\nThese are all the time values measured in hours. Are they correct? (Enter 1 for "Yes" or 0 for "No") ';
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      Pos = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	Pos = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    %% Find out where the subject is
    prompt = '\nWhich item in the list corresponds to the subject/condition? (format [1,2,3] or all) ';
    su_ind = input(prompt);
    [Subject_name,~,Subject] = unique({Info{su_ind,pr_ind}});

  else
    % New case where time points must be read from title field or somewhere else.      
    timePointsMatrix = ExtractTimePoints(dis(pr_ind));
      
    Pos = [];
    timeUnit = 'N/A';
      
    if not(isempty(timePointsMatrix))
      Pos = cell2mat(timePointsMatrix(:,1))';
      timeUnit = timePointsMatrix(1,2);
    end
	
    sane_check = 0;
    while sane_check == 0
      display(Pos');
      prompt = ['These are all the time values measured in ' timeUnit{1} '. Are they correct? (Enter 1 for "Yes" or 0 for "No") '];
      sane_check = input(prompt);
      if sane_check ~= 0
	break;
      end
      timePointsMatrix = InputTimePointsManually();
	
      Pos = [];
      timeUnit = 'N/A';
	
      if not(isempty(timePointsMatrix))
	Pos = cell2mat(timePointsMatrix(:,1))';
	timeUnit = timePointsMatrix(1,2);
      end
    end

    if iscellstr(str_ind)
      Subject_name = strjoin(str_ind, '_');
    else
      Subject_name = str_ind;
    end
    Subject = repmat(1, 1, size(Pos,2));
  end
    
  
