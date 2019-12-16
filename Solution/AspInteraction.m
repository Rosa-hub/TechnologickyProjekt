classdef AspInteraction
    properties
        
        propType
        streamName
        componentName
        xlsObj
        
    end
    
    properties(Constant)
        
        xlsfiel='\Connection\Aspeninteraction.xlsm'
        aspfile='\Connection\SPD.apwz'
        dataPath='\Data\Streams\'
        
    end
    
    methods
        
        function [obj,vargout] = openAspConnection(obj)
            
            try
                vargout=0;
                excelObject=actxserver('Excel.Application');
                obj.xlsObj=excelObject;
                excelObject.Visible = 0;
                excelObject.Workbooks.Open(fullfile(pwd,obj.xlsfiel));
                excelObject.Run('aspinit',fullfile(pwd,obj.aspfile));
                
                
                vargout=1;
                if vargout==0
                    error('AspInteraction:ConncectionFailed','Internal conncetion error')
                end
            catch e
                if vargout~=0
                e=MException('AspInteraction:ConncectionFailed','Unable to establish Aspen connection!');
                end
                throw(e)
            end
            
        end
        
        function vargout = closeAspConnection(obj)
            
            try
                vargout=0;
                obj.xlsObj.Run('aspClose');
                obj.xlsObj.DisplayAlerts = false;
                obj.xlsObj.Quit;
                obj.xlsObj.delete;
                
                
                vargout=1;
                if vargout==0
                    error('AspInteraction:ConncectionFailed','Internal conncetion error')
                end
                
            catch e
                
                if vargout~=0
                    e=MException('AspInteraction:ConncectionFailed','Unable to establish Aspen connection!');
                end
                
                throw(e)
                
            end
            
        end
        
        function aspRunSimulation(obj)
            obj.xlsObj.Run('aspRunSimulation')
        end
        
        function aspSetScalar(obj,val)
        
            fullpath=strcat(obj.dataPath,obj.streamName,'\Input\',obj.propType,'\MIXED\',obj.componentName);
            obj.xlsObj.Run('aspSetScalar',fullpath,val);
        end
        
        function vargout = aspGetScalar(obj)
            fullpath=strcat(obj.dataPath,obj.streamName,'\Output\',obj.propType,'\MIXED\',obj.componentName);
            vargout=obj.xlsObj.Run('aspGetScalar',fullpath);
        end
    end
    
end
    
    