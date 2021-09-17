function readPolars(this)
    this.Polars = {};
    for i=1:length(this.PolarFiles)
        try
            fid = fopen(this.PolarFiles{i});
            status = fseek(fid, 429, 'bof');
            data=textscan(fid,'%f%f%f%f%f%f%f');
            fclose(fid);
            if not(isempty(data)) && status == 0
                data=cell2mat(data);
            else
                data = [];
            end
            
            n=any(isnan(data),2);
            data(n,:) = [];
            Polar.Alpha   = data(:,1);
            Polar.CL      = data(:,2);
            Polar.CD      = data(:,3);
            Polar.CDp     = data(:,4);
            Polar.CM      = data(:,5);
            Polar.Top_Xtr = data(:,6);
            Polar.Bot_Xtr = data(:,7);
            
            this.Polars{i} = Polar;
        catch
            error('Invalid polar file %s', this.PolarFiles{i} );
        end
    end
    
    if ~this.KeepFiles 
        for i=1:length(this.PolarFiles)
            delete(this.PolarFiles{i});
        end
        this.PolarFiles = {};
        delete(this.ActionsFile);
        if isa(this.Airfoil,'Airfoil')
            delete(this.AirfoilFile);
        end
    end
    
end