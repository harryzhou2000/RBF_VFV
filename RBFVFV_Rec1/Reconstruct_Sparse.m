function elems = Reconstruct_Sparse(elems,A,Bw,options,levels)

if(nargin == 4)
    b = Bw*[elems.value]';
%     rec  =A\b;
    rec = F_CGSolve(A,b);
    ndof = elems(1).dofRec;
    for it = 1:numel(elems)
    	elems(it).recvalue(:) =  rec((it-1)*ndof+1:it*ndof);
    end
    
    
    
    %% Jacobi Iteration
    % while err > err0 && iter < itmax
    %     for it = 1:numel(elems)
    %         elems(it).recvalueN = elems(it).Ainv * ...
    %             (...
    %             elems(it).B{1}*elems(elems(it).neighbours(1)).recvalue...
    %             +elems(it).B{2}*elems(elems(it).neighbours(2)).recvalue...
    %             +elems(it).b...
    %             );
    %     end
    %     for it = 1:numel(elems)
    %         elems(it).recvalue = elems(it).recvalueN;
    %     end
    %     iter = iter+1;
    %
    %     if mod(iter,10) == 0
    %         err = 0;
    %         for it = 1:numel(elems)
    %             err = max(err,norm(elems(it).A*elems(it).recvalue - ...
    %                 elems(it).B{1}*elems(elems(it).neighbours(1)).recvalue -...
    %                 elems(it).B{2}*elems(elems(it).neighbours(2)).recvalue -...
    %                 elems(it).b , inf));
    %         end
    % %         fprintf('iter %d rec err = %e\n',iter,err);
    %     end
    % end
    % fprintf('End: %d rec err = %e\n',iter,err);
    % if iter == itmax
    %     warning('reconstruction not exact');
    % end
    
    %% Calculating Extented
    % variation
    %     for it = 1:numel(elems)
    %         elems(it).recvalueL = elems(it).AinvL*(...
    %             +elems(it).BL{1}*elems(elems(it).neighbours(1)).recvalue...
    %             +elems(it).BL{2}*elems(elems(it).neighbours(2)).recvalue...
    %             +elems(it).bL...
    %             );
    %
    %         elems(it).recvalueR = elems(it).AinvR*(...
    %             +elems(it).BR{1}*elems(elems(it).neighbours(1)).recvalue...
    %             +elems(it).BR{2}*elems(elems(it).neighbours(2)).recvalue...
    %             +elems(it).bR...
    %             );
    %     end
    switch options.extend
        case {1,3,4}
            %%%%%%%%% using lsq fitting 李万爱
            for it = 1:numel(elems)
                elems(it).recvalueL = elems(it).TlsqL * elems(it).recvalue;
                elems(it).recvalueR = elems(it).TlsqR * elems(it).recvalue;
            end
        case 2
            %%%%%%%%% using lsq fitting 李万爱++
            for it = 1:numel(elems)
                left = elems(it).neighbours(1);
                right = elems(it).neighbours(2);
                elems(it).recvalueL = elems(it).TlsqL * [elems(it).recvalue;elems(left).recvalue] ...
                    + elems(it).UlsqL * (elems(it).value-elems(left).value);
                elems(it).recvalueR = elems(it).TlsqR * [elems(it).recvalue;elems(right).recvalue] ...
                    + elems(it).UlsqR * (elems(it).value-elems(right).value);
            end
    end
else
    %% using WBAP
    b = Bw*[elems.value]';
%         rec  =A\b;
    rec = F_CGSolve(A,b);
    ndof = elems(1).dofRec;
    for it = 1:numel(elems)
    	elems(it).recvalue(:) =  rec((it-1)*ndof+1:it*ndof);
    end
    
    for ilv = numel(levels):-1:1
        lv = levels(ilv);
        if ilv ==1
            lvb = 1;
        else
            lvb = levels(ilv-1)+1;
        end
        %% Calculating Extented
        %%%%%%%%%%% using variation reconstruction extension method xx
        %         for it = 1:numel(elems)
        %             elems(it).recvalueL = elems(it).AinvL*(...
        %                 +elems(it).BL{1}*elems(elems(it).neighbours(1)).recvalue...
        %                 +elems(it).BL{2}*elems(elems(it).neighbours(2)).recvalue...
        %                 +elems(it).bL...
        %                 );
        %
        %             elems(it).recvalueR = elems(it).AinvR*(...
        %                 +elems(it).BR{1}*elems(elems(it).neighbours(1)).recvalue...
        %                 +elems(it).BR{2}*elems(elems(it).neighbours(2)).recvalue...
        %                 +elems(it).bR...
        %                 );
        %         end
        for repeat = lvb:lv
            if options.constbound
                elems(1).recvalue =  elems(1).recvalue*0;
                elems(end).recvalue =  elems(end).recvalue*0;
            end
            switch options.extend
                case {1,3,4}
                    %%%%%%%%% using lsq fitting 李万爱
                    for it = 1:numel(elems)
                        elems(it).recvalueL = elems(it).TlsqL * elems(it).recvalue;
                        elems(it).recvalueR = elems(it).TlsqR * elems(it).recvalue;
                    end
                case 2
                    %%%%%%%%% using lsq fitting 李万爱++
                    for it = 1:numel(elems)
                        left = elems(it).neighbours(1);
                        right = elems(it).neighbours(2);
                        elems(it).recvalueL = elems(it).TlsqL * [elems(it).recvalue;elems(left).recvalue] ...
                            + elems(it).UlsqL * (elems(it).value-elems(left).value);
                        elems(it).recvalueR = elems(it).TlsqR * [elems(it).recvalue;elems(right).recvalue] ...
                            + elems(it).UlsqR * (elems(it).value-elems(right).value);
                    end
            end
            
            %% Limiting
            
            for it = 1:numel(elems)
                left = elems(it).neighbours(1);
                right = elems(it).neighbours(2);
                elems(it).recvalue(lvb:lv,1) = ...
                    WBAP_L2([elems(it).recvalue(lvb:lv,1),elems(left).recvalueR(lvb:lv,1),elems(right).recvalueL(lvb:lv,1)],...
                    10);
            end
        end
    end
    
end
if options.constbound
    elems(1).recvalue =  elems(1).recvalue*0;
    elems(end).recvalue =  elems(end).recvalue*0;
end


end

