% function A = mycost(px1,py1,vx1,vy1,S1,ax1,ay1,len,px,py,S,v95,maxJump)
function [A,Av,Aa,As] = mycost(px1,py1,vx1,vy1,S1,ax1,ay1,len,px,py,S,v95,maxJump,countLost)
%MYCOST ����������Ϣ����ɱ�����
%   px1(t-1),py1(t-1),vx1(t-2),vy1(t-2),ax1(t-3),ay1(t-3),px(t),py(t)
%   Ҫ��[px(t),py(t)]��[px1(t-1),py1(t-1)]ƥ��
%   px1 M��1�У�MΪ������Ŀ�� px N��1�У� NΪʶ����ĸ�����Ŀ�� M��N���Բ����
%   ͨ�� N >= M, ��֡ʱ N < M.

%   2015/9/8  ���ڶ�ʧ�ĸ������³��ֵģ���������κ�λ���ǶԵģ�����������ɱ�����Ϊ0��ᵼ����ʶ�𵽵ĸ�������
%             ����ʵ�ʸ�������ʱ����ʧ�ĸ���ռ��δ��ʧ�����λ�ã���δ��ʧ�ĸ������µĶ�ʧ���塣
%             ��ˣ�A(ii,:)=9999999999(�㹻���������������Inf)��
%             �޸���Ծ������������֡��Խ�󣬿������������ΧԽ��
%             �����Ѿ���ʧ�ĸ���ҲҪ�޶�������Χ���ڿ������������Χ�ڵ�costһ����֮���Inf
%   2015/4/11 �޸�Ȩ�أ�������仯��Ȩ�ر�С��ֻ���ٶȺͼ��ٶȲ���Сʱ�ŷ������ã���Ϊ���������Ϊ��ָ���ɴ���

    if nargin<12
        v95=40;
    end
    if nargin<13
        maxJump=1;
    end
    A = zeros(size(px1,1),size(px,1));
    
                %%% test %%%
    Av = zeros(size(px1,1),size(px,1));
    Aa = zeros(size(px1,1),size(px,1));
    As = zeros(size(px1,1),size(px,1));
                %%% test %%%
                
    wvel = 10; % �ٶȳ��������ٶȵ�Ȩ�أ��ɴ��ų�Զ�ġ�
    wacce = 5; % ���ٶȵ�Ȩ�أ��ɴ��ų���������ģ���Ӧ�ÿ��Ǽ��ٶ���
    warea = 3; % �����Ȩ�أ�
    
%     vel = sqrt((vx1(:).^2 + vy1(:).^2)); % t-2ʱ���ٶȴ�С
    for ii = 1:size(px1,1)
        % �����δ��ʧ����ʧ���Ӷ�ʧ���ٴ��һ�ʱ��������ص�λ�ã��ٶȣ����ٶȶ�������NaN
        % ��˼���ʱ��Ҫ�������������ڵ���C���Գ���ʱ�ᵼ�³������
        if isnan(px1(ii)) || isnan(vx1(ii)) || isnan(ax1(ii))
            % �����������һ֡�Ѿ���ȫ��ʧ�ˣ������ڿ������κ�λ��
            A(ii,:) = 9999999999; % ������ɱ�ȫ��һ������Ϊ9999999999
            % 2015.9.8 �����Ѿ���ʧ�ĸ���ҲҪ�޶�������Χ���ڿ������������Χ�ڵ�costһ����֮���Inf
            for jj = 1:size(px,1)
                % �����ii��source���jj��targetƥ��ĳɱ����ɱ�ԽСԽ��
                vx2 = px(jj) - px1(ii); 
                vy2 = py(jj) - py1(ii);
                v2 = sqrt(vx2^2+vy2^2);  %ƥ���source��t-1ʱ���ٶȴ�С
                if v2 > v95 * (countLost(ii) + 2) + maxJump * len(ii) 
                    % 2015.9.8 �޸���Ծ������������֡��Խ�󣬿������������ΧԽ��
                    A(ii,jj) = A(ii,jj) + Inf;  % �������ٶ�̫�󣬱�����Ծ��
                end
            end
        else
            for jj = 1:size(px,1)
                % �����ii��source���jj��targetƥ��ĳɱ����ɱ�ԽСԽ��
                vx2 = px(jj) - px1(ii); 
                vy2 = py(jj) - py1(ii);
                v2 = sqrt(vx2^2+vy2^2);  %ƥ���source��t-1ʱ���ٶȴ�С
                ax2 = vx2 - vx1(ii); 
                ay2 = vy2 - vy1(ii);
                a2 = ax2^2 + ay2^2; % ƥ���source��t-2ʱ�̵ļ��ٶ�
                a1 = ax1(ii)^2 + ay1(ii)^2; % ƥ���source��t-3ʱ�̵ļ��ٶ�
                if v2 > v95
                    tmp = v2 - v95; % �ٶȳ���Խ�࣬�ɱ�Խ��
                else
                    tmp = 0;
                end
                
                %%% test %%%
                Av(ii,jj) = tmp;
                Aa(ii,jj) = sqrt(a1 + a2);
                As(ii,jj) = sqrt(abs(S(jj)-S1(ii)));
                %%% test %%%
                
                A(ii,jj) = wvel * tmp + wacce * sqrt(a1 + a2) + warea * sqrt(abs(S(jj)-S1(ii)));
                if v2 > v95 * (countLost(ii) + 2) + maxJump * len(ii) 
                    % 2015.9.8 �޸���Ծ������������֡��Խ�󣬿������������ΧԽ��
                    A(ii,jj) = A(ii,jj) + Inf;  % �������ٶ�̫�󣬱�����Ծ��
                end
            end
        end
    end
end