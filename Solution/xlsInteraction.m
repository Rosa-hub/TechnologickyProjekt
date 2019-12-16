classdef xlsInteraction
    
    properties
        file
        xlsObject
        WB
        Alsht
        Sht
    end
    
    properties(Constant)
        
        Dirrectory='Connection\'
        
    end
    
    methods
        
        function obj =xlsOpenConnection(obj,sheet)
            
            try
                filepath=strcat(obj.Dirrectory,obj.file);
                excelObject=actxserver('Excel.Application');
                excelObject.Visible = 0;
                obj.xlsObject.DisplayAlerts = false;

                obj.xlsObject=excelObject;
                obj.WB=excelObject.Workbooks.Open(fullfile(pwd,filepath));
                obj.Alsht=obj.WB.Sheets;

                if exist('sheet','var') ==1
                    obj.Sht = get(obj.Alsht,'Item',sheet);
                else
                    obj.Sht = get(obj.Alsht,'Item',obj.Alsht.Count);    
                end

            catch e
                obj=obj.xlsCloseConnection;
                rethrow(e)
            end
            
            
                
        end
        
        function EconomyInitialize(obj)
            x=xlsInteraction;
            x.file='Economy.xlsm';
            x=x.xlsOpenConnection;
            x.xlsObject.Visible=1;
            x.xlsObject.Run('copyTemplate');   
        end
        
        function xlsRunMacro(obj,macro,params)
            try
                if exist('params','var')==0
                    obj.xlsObject.Run(macro)
                else
                    obj.xlsObject.Run(macro,params)
                end
            catch e
                obj=obj.xlsCloseConnection;
                rethrow(e)
            end
        end
        
        function obj=xlsCloseConnection(obj)
            
            obj.xlsObject.DisplayAlerts = false;
            obj.WB.Save
            obj.xlsObject.Quit;
            obj.xlsObject.delete;
         
            obj=xlsInteraction;
            
        end
        
        function xlsDataWrite(obj,range,val)
           try 
                RngObj=obj.Sht.Range(range);
                RngObj.Value = val;
            
           catch e
               
               obj=obj.xlsCloseConnection;
               rethrow(e)
           end
                        
        end
        
        function out = xlsDataRead(obj,range)
           try
            
               RngObj=obj.Sht.Range(range); 
               out=RngObj.Value;
               
               if iscell(out)==1
                   out=cell2mat(out);
               end
               
           catch e
               obj=obj.xlsCloseConnection;
               rethrow(e)
           end
        end
    end
    
end
