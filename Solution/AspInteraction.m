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
        
        function [obj,vargout] = openConnectionExcel(obj)
            
            try
                vargout=0;
                excelObject=actxserver('Excel.Application');
                obj.xlsObj=excelObject;
                excelObject.Visible = 1;                
                
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
        
        function obj = openAspConnection(obj)
            obj=obj.openConnectionExcel;            
            obj.xlsObj.Workbooks.Open(fullfile(pwd,obj.xlsfiel));
            obj.xlsObj.Run('openConnection');
        end
        
        function vargout = aspCloseAspConnection(obj)
            
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
        
        function aspSetScalar(obj,val,custpath)
            if isempty(obj.componentName)
                hash='';
            else
                hash='\';
            end
            if exist('custpath','var')==0
            fullpath=strcat(obj.dataPath,obj.streamName,'\Input\',obj.propType,'\MIXED',hash,obj.componentName);
            else
              fullpath=strcat(custpath);  
            end
            obj.xlsObj.Run('aspSetScalar',fullpath,val);
        end
        
        function vargout = aspGetScalar(obj,custpath)
            if isempty(obj.componentName)
                hash='';
            else
                hash='\';
            end
                
            if exist('custpath','var')==0
               fullpath=strcat(obj.dataPath,obj.streamName,'\Output\',obj.propType,'\MIXED',hash,obj.componentName);
            else
              fullpath=strcat(custpath);  
            end
            vargout=obj.xlsObj.Run('aspGetScalar',fullpath);
        end
        
        function obj=clearParams(obj)
        obj.propType=[];
        obj.streamName=[];
        obj.componentName = [];
        end
        
    end
    

    
end
    
    