classdef xlsInteraction
    
    properties
        file
        xlsObject
        WB
        Alsht
        sheet

    end
    
    properties(Constant)
        
        Dirrectory='Connection\'
        
    end
    
    methods
        
        function obj =xlsOpenConnection(obj)
            
            try
                filepath=strcat(obj.Dirrectory,obj.file);
                excelObject=actxserver('Excel.Application');
                excelObject.Visible = 0;
                excelObject.DisplayAlerts = false;

                obj.xlsObject=excelObject;
                obj.WB=excelObject.Workbooks.Open(fullfile(pwd,filepath));
                obj.Alsht=obj.WB.Sheets;

            catch e
                obj=obj.xlsCloseConnection;
                rethrow(e)
            end
            
            
                
        end
        
        function EconomyInitialize(obj)
            x=xlsInteraction;
            x.file='Economy_test.xlsm';
            x=x.xlsOpenConnection;
            x.xlsObject.Visible=1;
            x.xlsRunMacro('copyTemplate');  
            x.xlsDataWrite('A5:A10',20);
            x=x.xlsCloseConnection;
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
            obj.WB.Save;
            obj.xlsObject.Quit;
            obj.xlsObject.delete;
            
        end
        
        function xlsDataWrite(obj,range,val)
           try 
               
                if isempty(obj.sheet)==0
                    Sht = get(obj.Alsht,'Item',obj.sheet);
                else
                    Sht = get(obj.Alsht,'Item',obj.Alsht.Count);    
                end
               
                RngObj=Sht.Range(range);
                RngObj.Value = val;
            
           catch e
               
               obj=obj.xlsCloseConnection;
               rethrow(e)
           end
                        
        end
        
        function out = xlsDataRead(obj,range)
           try
               
                if isempty(obj.sheet) ==0
                    Sht = get(obj.Alsht,'Item',obj.sheet);
                else
                    Sht = get(obj.Alsht,'Item',obj.Alsht.Count);    
                end
            
               RngObj=Sht.Range(range); 
               out=RngObj.Value;
               
               if iscell(out)==1
                   out=cell2mat(out);
               end
               
           catch e
               obj=obj.xlsCloseConnection;
               rethrow(e)
           end
        end
        
        function out = getRange(obj,ref,val)
            
           
           Sht = get(obj.Alsht,'Item',obj.Alsht.Count);
           Rng1 = Sht.Range(ref);
           [a,b]=size(val);
           ref1=Rng1.Address;
           ref2=Rng1.get('Offset',a-1,b-1).Address;
           
           out=strcat(ref1,":",ref2);
            
        end
    end
    
end
