clc
close all

angles = [];
headX = [];
headY = [];
X2 = [];
Y2 = [];
X3 = [];
Y3 = [];
centroidX = [];
centroidY = [];
X5 = [];
Y5 =[];
X6 = [];
Y6 =[];
tailX = [];
tailY = [];
TSdis = []; % �ڶ����ο���͵������ο���֮��ľ���
htdistance = [];
hDistances = [];
vector_shX = [];
vector_shY = [];
vector_ccX = [];
vector_ccY = [];
vector_hhX = [];
vector_hhY = [];
vector_ttX = [];
vector_ttY = [];
sh_hh_dot_product = [];
sh_tt_dot_product = [];
sh_cc_dot_product = [];
cc_hh_dot_product = [];

trainPath = 'D:\Postgraduate\NewDL\data\N2\5\Binary_85\';
% trainPath = 'D:\Postgraduate\NewDL\escape\19\Binary\';
theFiles = dir([trainPath, '*.jpg']);
disp(length(theFiles));
count = 1;
train_num = length(theFiles);

sort_nat_name = sort_nat({theFiles.name});   % �������ݼ������ݵ��������� �����ݽ�������

for k = 1:train_num
    fullFileName = sort_nat_name{k};
    fprintf(1, 'Now reading %s\n', fullFileName);
    I = imread([trainPath, fullFileName]);
    
    if k == 1
        [height, width] = size(I);
        % ��ȡͼ�������
        [~, imgName, ~] = fileparts(fullFileName);
        [headx,heady,x2,y2,x3,y3,centroidx,centroidy,x5,y5,x6,y6,tailx,taily,TSdistancs] = plotDivideSpline3(I, imgName);
        figure('Name', '9_Skeleton');
        
        angle = headBendAngle1(headx,heady,x2,y2,x3,y3);
        fprintf('head bend angle is %.2f\n', angle);
        angles = [angles, angle];
        
        heady = height - heady;  % ��ԭ
        taily = height - taily;  
        headToTailDist = ptDist(headx, heady, tailx, taily);
        hdistance = point_to_line_segment_distance(centroidx, centroidy, headx, heady, tailx, taily);
        
        hDistances = [hDistances, hdistance];
        htdistance = [htdistance, headToTailDist];
        TSdis = [TSdis, TSdistancs];
        headX = [headX,headx];
        headY = [headY,heady];
        X2 = [X2, x2];
        Y2 = [Y2, y2];
        X3 = [X3, x3];
        Y3 = [Y3, y3];
        centroidX = [centroidX, centroidx];
        centroidY = [centroidY, centroidy]; 
        X5 = [X5, x5];
        Y5 = [Y5, y5];
        X6 = [X6, x6];
        Y6 = [Y6, y6];
        tailX = [tailX, tailx];
        tailY = [tailY, taily];
        vector_shX = [vector_shX, headx-centroidx];
        vector_shY = [vector_shY, heady-centroidy];
        vector_ccX = [vector_ccX, 0];
        vector_ccY = [vector_ccY, 0];
        vector_hhX = [vector_hhX, 0];
        vector_hhY = [vector_hhY, 0];
        vector_ttX = [vector_ttX, 0];
        vector_ttY = [vector_ttY, 0];
        sh_hh_dot_product = [sh_hh_dot_product, 0];
        sh_tt_dot_product = [sh_tt_dot_product, 0];
        sh_cc_dot_product = [sh_cc_dot_product, 0];
        cc_hh_dot_product = [cc_hh_dot_product, 0];
        continue;
    end
    
    [newHeadx,newHeady,x2,y2,x3,y3,centroidx,centroidy,x5,y5,x6,y6,newTailx,newTaily,ang,dis,TSdistancs,hdistance,hhx,hhy,ttx,tty] = findNextHead3(I, headx, heady, tailx, taily);
    
    % ������ر���
    angles = [angles, ang];
    headx = newHeadx;
    heady = newHeady;
    tailx = newTailx;
    taily = newTaily;
    
    headX = [headX,headx];
    headY = [headY,heady];
    X2 = [X2, x2];
    Y2 = [Y2, y2];
    X3 = [X3, x3];
    Y3 = [Y3, y3];
    centroidX = [centroidX, centroidx];
    centroidY = [centroidY, centroidy]; 
    X5 = [X5, x5];
    Y5 = [Y5, y5];
    X6 = [X6, x6];
    Y6 = [Y6, y6];
    tailX = [tailX, tailx];
    tailY = [tailY, taily];
    htdistance = [htdistance, dis];
    TSdis = [TSdis, TSdistancs];
    hDistances = [hDistances, hdistance];
    vector_shX = [vector_shX, headx-centroidx];
    vector_shY = [vector_shY, heady-centroidy];
    vector_hhX = [vector_hhX, hhx];
    vector_hhY = [vector_hhY, hhy];
    vector_ttX = [vector_ttX, ttx];
    vector_ttY = [vector_ttY, tty];
end

% �쳣�ж�
noChangeCount = 0;  % ��¼����û�иı��ͼƬ����
exceptionIndex = 0;  % ��¼�쳣���ֵ�ͼƬ����

% ����ͷ���ڶ������
anglength = length(angles);
newangles = [];
i = 1;

% ѭ���������½Ƕ�����newangles�н�ֻ��������ֵ����0�ĽǶȡ�
for t = 1:(anglength - 1)
    if abs(angles(t)) > 5
        newangles(i) = angles(t);
        i = i + 1;
    end
end

a = length(newangles);
anglecha = [];

% ����������֮֡��ǶȵĲ�ֵ������anglecha
for t = 1:(a - 1)
    anglecha(t) = newangles(t) - newangles(t + 1);
end
% ���Ե���������֮֡��İڶ��ǶȲ�ֵС��5��Ľǣ�����������ĽǶȵ�λ��
x = find(abs(anglecha) < 5);
anglecha(x) = [];
headThrashesnum = 0;

headThrashes = {};  % �洢����ͷ��ҡ�ε�ͼƬ��
for i = 1:(length(anglecha) - 1)
    if anglecha(i) > 0 && anglecha(i + 1) < 0
        if abs(abs(anglecha(i)) - abs(anglecha(i + 1))) > 5
            headThrashesnum = headThrashesnum + 1;
            imgName = sort_nat_name{i + 1};  % ��ȡ����ͷ��ҡ�ε�ͼƬ��
            headThrashes{end+1} = imgName;  % ��ͼƬ����ӵ��б���
        end
    end
end
% �����ص��������� ����һ�ΰڶ�
headThrashesnum = round(headThrashesnum / 2);

% �������ͷ��ҡ�ε�ͼƬ��
fprintf('Image names with head thrashes:\n');
for i = 1:length(headThrashes)
    fprintf('%s\n', headThrashes{i});
end
fprintf('The number of head thrash are %d\n', headThrashesnum)

for k = 1:train_num-1
    headToTailDist1 = htdistance(k);  % ��ȡ��ǰѭ���е�ͷβ����
    headToTailDist2 = htdistance(k+1);  % ��ȡ��ǰѭ���е�ͷβ����
    
    % �ж��Ƿ�����omega turn������
    if noChangeCount < 50
        if headToTailDist1 == headToTailDist2
            noChangeCount = noChangeCount + 1;
        else
            noChangeCount = 0;
        end 
    else
        % disp(['�쳣�����ڣ�', sort_nat_name{k-50}]);
    end
end

if noChangeCount < 50
    disp('û�з����쳣');
end

omegaTurns = 0;
omegaTurnStart = false;  % omega turn��ʼ���
OmegaTurnFrame = [];
maxHeadToTailDist = max(htdistance);  % ���ͷβ����
fprintf('maxHeadToTailDist is %.2f\n', maxHeadToTailDist);
omegaTurnImgs = {};  % ���ڱ�������omega turn������ͼ������

for k = 1:train_num
    headToTailDist = htdistance(k);  % ��ȡ��ǰѭ���е�ͷβ����
    % �ж��Ƿ�����omega turn������
    if (headToTailDist < maxHeadToTailDist * 0.5) || (headToTailDist < maxHeadToTailDist * 0.6  && TSdis(k) < maxHeadToTailDist * 0.4)
        omegaTurnStart = true;
        disp(['����������ͼ�����ƣ�', sort_nat_name{k}, '��headToTailDistֵ��', num2str(headToTailDist)]);
        OmegaTurnFrame = [OmegaTurnFrame,k];
    else 
        if omegaTurnStart
            omegaTurns = omegaTurns + 1;  % ����omega turn����
            disp(['������ͼ�����ƣ�', sort_nat_name{k}, '��headToTailDistֵ��', num2str(headToTailDist)]);
            omegaTurnStart = false;
        end
    end
end

% �����������
disp([OmegaTurnFrame]);
fprintf('The number of omega turns are %d\n', omegaTurns)

for i = 2:train_num
    vector_ccX(i) = centroidX(i) - centroidX(i-1);
    vector_ccY(i) = centroidY(i) - centroidY(i-1);
end

% ѭ��������
for i = 2:train_num
    % ���� sh ������ hh ֮��ĵ��
    sh_hh_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_hhX(i); vector_hhY(i)]);
    
    % ���� sh ������ tt ֮��ĵ��
    sh_tt_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_ttX(i); vector_ttY(i)]);
    
    % ���� sh ������ tt ֮��ĵ��
    sh_cc_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_ccX(i); vector_ccY(i)]);
    
    % ���� cc ������ hh ֮��ĵ��
    cc_hh_dot_product(i) = dot([vector_ccX(i); vector_ccY(i)], [vector_hhX(i); vector_hhY(i)]);
    
end

reversalNum = 0;
reversals = [];
reversalFrame = []; % �洢����������֡������

% ѭ��������������Ƿ�С����
consecutive_count = 0; % ��¼����С�������ݵĸ���
for i = 1:train_num
    if sh_hh_dot_product(i) < 0 && sh_tt_dot_product(i) < 0 && sh_cc_dot_product(i) <0 
        consecutive_count = consecutive_count + 1;
    else
        % �����������������������
        consecutive_count = 0;
    end
    reversals(i) = consecutive_count;
end

for i = 1:train_num-1
    % �������С�������ݸ����ﵽ3�������ϣ��������1������¼֡����
    if reversals(i) >= 3 && reversals(i+1)==0
        reversalNum = reversalNum + 1;
        reversalFrame(reversalNum) =  i ; % ��¼���һ������������֡������
    end
end

%{
for i = 1:reversalNum-1
    % �����10֡�ڱ�ǹ�һ��reversal���ͺ��Դ˴α��
    if reversalFrame(i+1) < (reversalFrame(i) + 5)
        reversalNum = reversalNum - 1;
    end
end
%}

% �����������
fprintf('The number of reversals are %d\n', reversalNum)
% �������������֡������
disp(reversalFrame);

PirouetteFrames = []; % ��ʼ���洢����������֡������
count = 0; % ��ʼ��������

for i = 1:length(reversalFrame)
    % ����OmegaTurnFrame���뵱ǰreversalFrame���ڵ�֡
    nearbyOmegaTurnFrames = OmegaTurnFrame(OmegaTurnFrame >= reversalFrame(i) - 15 & OmegaTurnFrame <= reversalFrame(i) + 15);
    
    if ~isempty(nearbyOmegaTurnFrames)
        % ������������֡��ӵ�PirouetteFrames�У��������ظ���֡
        PirouetteFrames = [PirouetteFrames, nearbyOmegaTurnFrames];
    end
end

% ����ȥ���ظ�֡
PirouetteFrames = unique(PirouetteFrames);

for i = 1:length(PirouetteFrames)-1
    % �����ǰ֡����һ֡�������������Ӽ�����
    if PirouetteFrames(i) + 1 ~= PirouetteFrames(i+1)
        count = count + 1;  
    end
end

if PirouetteFrames ~= 0
    count = count +1;
end

disp(['����Pirouette������֡��', num2str(PirouetteFrames)]);
disp(['The number of pirouettes are ', num2str(count)]);


% �������
data = table(sort_nat_name',headX',headY',X2',Y2',X3',Y3',centroidX',centroidY',X5',Y5',X6',Y6',tailX',tailY',htdistance',hDistances',sh_hh_dot_product',sh_tt_dot_product',sh_cc_dot_product', 'VariableNames', {'ImageName','headX','headY','X2','Y2','X3','Y3','centroidX','centroidY','X5','Y5','X6','Y6','tailX','tailY', 'HeadToTailDistance','hDistances','sh_hh_dot_product','sh_tt_dot_product','sh_cc_dot_product'});

% ȡ������
data.centroidX = round(data.centroidX);
data.centroidY = round(data.centroidY);

% �����д��CSV�ļ�
csvFileName = 'D:\Postgraduate\NewDL\data\N2\5\data\new_data.csv';
writetable(data, csvFileName);

