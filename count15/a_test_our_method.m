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
TSdis = []; % 第二个参考点和第六个参考点之间的距离
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

sort_nat_name = sort_nat({theFiles.name});   % 按照数据集中数据的命名规律 对数据进行排序

for k = 1:train_num
    fullFileName = sort_nat_name{k};
    fprintf(1, 'Now reading %s\n', fullFileName);
    I = imread([trainPath, fullFileName]);
    
    if k == 1
        [height, width] = size(I);
        % 获取图像的名称
        [~, imgName, ~] = fileparts(fullFileName);
        [headx,heady,x2,y2,x3,y3,centroidx,centroidy,x5,y5,x6,y6,tailx,taily,TSdistancs] = plotDivideSpline3(I, imgName);
        figure('Name', '9_Skeleton');
        
        angle = headBendAngle1(headx,heady,x2,y2,x3,y3);
        fprintf('head bend angle is %.2f\n', angle);
        angles = [angles, angle];
        
        heady = height - heady;  % 还原
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
    
    % 更新相关变量
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

% 异常判断
noChangeCount = 0;  % 记录连续没有改变的图片数量
exceptionIndex = 0;  % 记录异常出现的图片索引

% 计算头部摆动的情况
anglength = length(angles);
newangles = [];
i = 1;

% 循环结束后，新角度序列newangles中将只包含绝对值大于0的角度。
for t = 1:(anglength - 1)
    if abs(angles(t)) > 5
        newangles(i) = angles(t);
        i = i + 1;
    end
end

a = length(newangles);
anglecha = [];

% 将连续的两帧之间角度的差值保存在anglecha
for t = 1:(a - 1)
    anglecha(t) = newangles(t) - newangles(t + 1);
end
% 忽略掉连续的两帧之间的摆动角度差值小于5°的角，返回有意义的角度的位置
x = find(abs(anglecha) < 5);
anglecha(x) = [];
headThrashesnum = 0;

headThrashes = {};  % 存储发生头部摇晃的图片名
for i = 1:(length(anglecha) - 1)
    if anglecha(i) > 0 && anglecha(i + 1) < 0
        if abs(abs(anglecha(i)) - abs(anglecha(i + 1))) > 5
            headThrashesnum = headThrashesnum + 1;
            imgName = sort_nat_name{i + 1};  % 获取发生头部摇晃的图片名
            headThrashes{end+1} = imgName;  % 将图片名添加到列表中
        end
    end
end
% 将来回的弯曲动作 记作一次摆动
headThrashesnum = round(headThrashesnum / 2);

% 输出发生头部摇晃的图片名
fprintf('Image names with head thrashes:\n');
for i = 1:length(headThrashes)
    fprintf('%s\n', headThrashes{i});
end
fprintf('The number of head thrash are %d\n', headThrashesnum)

for k = 1:train_num-1
    headToTailDist1 = htdistance(k);  % 获取当前循环中的头尾距离
    headToTailDist2 = htdistance(k+1);  % 获取当前循环中的头尾距离
    
    % 判断是否满足omega turn的条件
    if noChangeCount < 50
        if headToTailDist1 == headToTailDist2
            noChangeCount = noChangeCount + 1;
        else
            noChangeCount = 0;
        end 
    else
        % disp(['异常发生在：', sort_nat_name{k-50}]);
    end
end

if noChangeCount < 50
    disp('没有发生异常');
end

omegaTurns = 0;
omegaTurnStart = false;  % omega turn开始标记
OmegaTurnFrame = [];
maxHeadToTailDist = max(htdistance);  % 最大头尾距离
fprintf('maxHeadToTailDist is %.2f\n', maxHeadToTailDist);
omegaTurnImgs = {};  % 用于保存满足omega turn条件的图像名称

for k = 1:train_num
    headToTailDist = htdistance(k);  % 获取当前循环中的头尾距离
    % 判断是否满足omega turn的条件
    if (headToTailDist < maxHeadToTailDist * 0.5) || (headToTailDist < maxHeadToTailDist * 0.6  && TSdis(k) < maxHeadToTailDist * 0.4)
        omegaTurnStart = true;
        disp(['满足条件的图像名称：', sort_nat_name{k}, '，headToTailDist值：', num2str(headToTailDist)]);
        OmegaTurnFrame = [OmegaTurnFrame,k];
    else 
        if omegaTurnStart
            omegaTurns = omegaTurns + 1;  % 增加omega turn计数
            disp(['结束的图像名称：', sort_nat_name{k}, '，headToTailDist值：', num2str(headToTailDist)]);
            omegaTurnStart = false;
        end
    end
end

% 输出计数次数
disp([OmegaTurnFrame]);
fprintf('The number of omega turns are %d\n', omegaTurns)

for i = 2:train_num
    vector_ccX(i) = centroidX(i) - centroidX(i-1);
    vector_ccY(i) = centroidY(i) - centroidY(i-1);
end

% 循环计算点积
for i = 2:train_num
    % 向量 sh 和向量 hh 之间的点积
    sh_hh_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_hhX(i); vector_hhY(i)]);
    
    % 向量 sh 和向量 tt 之间的点积
    sh_tt_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_ttX(i); vector_ttY(i)]);
    
    % 向量 sh 和向量 tt 之间的点积
    sh_cc_dot_product(i) = dot([vector_shX(i); vector_shY(i)], [vector_ccX(i); vector_ccY(i)]);
    
    % 向量 cc 和向量 hh 之间的点积
    cc_hh_dot_product(i) = dot([vector_ccX(i); vector_ccY(i)], [vector_hhX(i); vector_hhY(i)]);
    
end

reversalNum = 0;
reversals = [];
reversalFrame = []; % 存储满足条件的帧的索引

% 循环检查连续数据是否小于零
consecutive_count = 0; % 记录连续小于零数据的个数
for i = 1:train_num
    if sh_hh_dot_product(i) < 0 && sh_tt_dot_product(i) < 0 && sh_cc_dot_product(i) <0 
        consecutive_count = consecutive_count + 1;
    else
        % 如果不连续，重置连续计数
        consecutive_count = 0;
    end
    reversals(i) = consecutive_count;
end

for i = 1:train_num-1
    % 如果连续小于零数据个数达到3个及以上，则计数加1，并记录帧索引
    if reversals(i) >= 3 && reversals(i+1)==0
        reversalNum = reversalNum + 1;
        reversalFrame(reversalNum) =  i ; % 记录最后一个满足条件的帧的索引
    end
end

%{
for i = 1:reversalNum-1
    % 如果在10帧内标记过一次reversal，就忽略此次标记
    if reversalFrame(i+1) < (reversalFrame(i) + 5)
        reversalNum = reversalNum - 1;
    end
end
%}

% 输出计数次数
fprintf('The number of reversals are %d\n', reversalNum)
% 输出满足条件的帧的索引
disp(reversalFrame);

PirouetteFrames = []; % 初始化存储满足条件的帧的数组
count = 0; % 初始化计数器

for i = 1:length(reversalFrame)
    % 查找OmegaTurnFrame中与当前reversalFrame相邻的帧
    nearbyOmegaTurnFrames = OmegaTurnFrame(OmegaTurnFrame >= reversalFrame(i) - 15 & OmegaTurnFrame <= reversalFrame(i) + 15);
    
    if ~isempty(nearbyOmegaTurnFrames)
        % 将满足条件的帧添加到PirouetteFrames中，并包括重复的帧
        PirouetteFrames = [PirouetteFrames, nearbyOmegaTurnFrames];
    end
end

% 排序并去除重复帧
PirouetteFrames = unique(PirouetteFrames);

for i = 1:length(PirouetteFrames)-1
    % 如果当前帧与下一帧不连续，则增加计数器
    if PirouetteFrames(i) + 1 ~= PirouetteFrames(i+1)
        count = count + 1;  
    end
end

if PirouetteFrames ~= 0
    count = count +1;
end

disp(['满足Pirouette条件的帧：', num2str(PirouetteFrames)]);
disp(['The number of pirouettes are ', num2str(count)]);


% 创建表格
data = table(sort_nat_name',headX',headY',X2',Y2',X3',Y3',centroidX',centroidY',X5',Y5',X6',Y6',tailX',tailY',htdistance',hDistances',sh_hh_dot_product',sh_tt_dot_product',sh_cc_dot_product', 'VariableNames', {'ImageName','headX','headY','X2','Y2','X3','Y3','centroidX','centroidY','X5','Y5','X6','Y6','tailX','tailY', 'HeadToTailDistance','hDistances','sh_hh_dot_product','sh_tt_dot_product','sh_cc_dot_product'});

% 取整操作
data.centroidX = round(data.centroidX);
data.centroidY = round(data.centroidY);

% 将表格写入CSV文件
csvFileName = 'D:\Postgraduate\NewDL\data\N2\5\data\new_data.csv';
writetable(data, csvFileName);

