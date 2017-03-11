#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>
#include <queue>
#include <set>
#include <cmath>

using namespace std;

int targetCol = 6; // 目标列号
vector<int> l_feature; // feature列表
set<int> l_used; //存放在Tree中已分割的feature
map<int, int> d_kv_fd; // 函数依赖对{A:B}
map<int, int> d_kv_fd_2; // 反序函数依赖对{B:A}
// 函数依赖对中后置值对应的数量{A:{Ai:{Bi:numb, Bi:numb}, Ai:{Bi:numb, Bi:numb}}, A:{Ai:{Bi:numb}}}
map<int, map<int, map<int, int> > > d_kv_fd_numb;
map<int, vector<int> > d_feature_numb; // 每个特征值中所包含的数值列表{A:[Ai, Ai], B:[Bi, Bi]}
map<int, map<int, int> > d_maxBi; // 存储Ai下最大的Bi值{A:{Ai:Bi, Ai:Bi}, A:{Ai:Bi, Ai:Bi}}

int data_max_size = 0;
vector<vector<int> > train_data, test_data;

double random(double start, double end){
    return start+(end-start)*rand()/(RAND_MAX + 1.0);
}

//读取数据至vector<vector<int>>中 adrr中的特征值，需要是离散字符变量转化为离散数值变量后的结果
vector<vector<int> > read_data(string adrr){

	//  l_feature  d_kv_fd  d_kv_fd_2
	for(int i=0; i<targetCol; i++)
		l_feature.push_back(i);
	///////////////////////////////////////////////////////////////////////////
//	d_kv_fd[1] = 2;
//	d_kv_fd_2[2] = 1;
//  d_kv_fd[3] = 4;
//	d_kv_fd_2[4] = 3;

    d_kv_fd[2] = 3;
    d_kv_fd_2[3] = 2;
	///////////////////////////////////////////////////////////////////////////

	// 读取数据
	ifstream file(adrr.c_str());
	string line;
	vector<vector<int> > data;
	while (getline(file, line)){   //整行读取，换行符“\n”区分，遇到文件尾标志eof终止读取
        istringstream sin(line); //将整行字符串line读入到字符串流istringstream中
        vector<int> fields; //声明一个字符串向量
        string field;
        while (getline(sin, field, ',')){ //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符
            fields.push_back(atoi(field.c_str())); //将刚刚读取的字符串添加到向量fields中
        }
        if(random(0, 9) > 5 && data.size() <= data_max_size)
            data.push_back(fields);
    }
    return data;
}

void pre_data(vector<vector<int> > data){
	long long dataSize = data.size();

    srand(unsigned(time(0)));

    //抽取1/3作为测试数据
    for(long long i=0; i<dataSize; i++){
    	if(int(random(0,9)) > 5)
    		test_data.push_back(data[i]);
    	else
    		train_data.push_back(data[i]);
    }

    long long trainDataSize = train_data.size();
    for(long long iData=0; iData<trainDataSize; iData++){
    	for(int iFeature=0; iFeature<l_feature.size(); iFeature++){


    		//d_kv_fd_numb构造
			map<int, int>::iterator kvIterator = d_kv_fd.find(l_feature[iFeature]);

			// l_feature[iFeature]是前置值A
    		if(kvIterator != d_kv_fd.end()){

				int A = l_feature[iFeature];
    			int Ai = train_data[iData][l_feature[iFeature]];
    			int Bi = train_data[iData][d_kv_fd[l_feature[iFeature]]];

    			map<int, map<int, map<int, int> > >::iterator AIterator = d_kv_fd_numb.find(A);
    			if(AIterator != d_kv_fd_numb.end()){

    				map<int, map<int, int> >::iterator AiIterator = d_kv_fd_numb[A].find(Ai);
    				if(AiIterator != d_kv_fd_numb[A].end()){

    					map<int, int>::iterator BiIterator = d_kv_fd_numb[A][Ai].find(Bi);
    					if(BiIterator != d_kv_fd_numb[A][Ai].end()){
    						int numb = d_kv_fd_numb[A][Ai][Bi];
    						d_kv_fd_numb[A][Ai][Bi] = numb+1;
    					}
    					else{
    						d_kv_fd_numb[A][Ai][Bi] = 1;
    					}

    				}
    				else{
    					map<int, int> BiTemp;
    					BiTemp[Bi] = 1;
    					d_kv_fd_numb[A][Ai] = BiTemp;
    				}
    			}
    			else{
    				map<int, int> BiTemp;
    				BiTemp[Bi] = 1;
    				map<int, map<int, int> > AiTemp;
    				AiTemp[Ai] = BiTemp;
    				d_kv_fd_numb[A] = AiTemp;
    			}
    		}

    		//d_feature_numb构造
    		map<int, vector<int> >::iterator featureIterator = d_feature_numb.find(l_feature[iFeature]);
    		if(featureIterator != d_feature_numb.end()){
    			int A = l_feature[iFeature];
    			int Ai = train_data[iData][A];
    			vector<int>::iterator AiIterator = find(d_feature_numb[A].begin(),d_feature_numb[A].end(),Ai);
    			if(AiIterator == d_feature_numb[A].end()){
    				d_feature_numb[A].push_back(Ai);
    			}
    		}
    		else{
    			int Ai = train_data[iData][l_feature[iFeature]];
    			vector<int> AiTemp;
    			AiTemp.push_back(Ai);
    			d_feature_numb[l_feature[iFeature]] = AiTemp;
    		}
    	}
    }

	for(map<int, map<int, map<int, int> > >::iterator AIterator=d_kv_fd_numb.begin(); AIterator!=d_kv_fd_numb.end(); ++AIterator){

		int A = AIterator->first;
		map<int, map<int, int> > AiMap = AIterator->second;

		for(map<int, map<int, int> >::iterator AiIerator=AiMap.begin(); AiIerator!=AiMap.end(); AiIerator++){
			int Ai = AiIerator->first;
			map<int, int> BiMap = AiIerator->second;

			int maxNumb = 0, maxBi = 0;

			for(map<int, int>::iterator BiIerator=BiMap.begin(); BiIerator!=BiMap.end(); BiIerator++){
				int Bi = BiIerator->first;
				int BiNumb = BiIerator->second;
				if(BiNumb > maxNumb){
					maxNumb = BiNumb;
					maxBi = Bi;
				}
			}
			d_maxBi[A][Ai] = maxBi;
		}
	}
}

const long long maxNodeNumb = 10000000;
long long treeNodeIndex = 0;
struct Node{
	map<int, Node*> subNode;
	vector<vector<int>* > nodeData;
	int splitNum;
	int maxClass;
	int treeLevel;
    long long unique;

	Node(){
		splitNum = -1;
		maxClass = -1;
		treeLevel = 0;
	}
    Node(int splitNumTemp, int maxClassTemp, int treeLevelTemp, long long uniqueTemp){
		splitNum = splitNumTemp;
		maxClass = maxClassTemp;
		treeLevel = treeLevelTemp;
        unique = uniqueTemp;
	}
}treeNode[maxNodeNumb];

int check_node_numb(Node* root){
    Node* tempNode = root;
    queue<Node*> queueNode;
    queueNode.push(tempNode);
    int nodeNumb = 0;
    while(!queueNode.empty()){
        tempNode = queueNode.front();
        queueNode.pop();
        nodeNumb++;
        for(map<int, Node*>::iterator it=tempNode->subNode.begin(); it!=tempNode->subNode.end(); it++)
            queueNode.push(it->second);
    }
    return nodeNumb;
}

//设置节点类所占最多值
int set_class(vector<vector<int>* > data){
	map<int, int> ans;
	for(int i=0; i<data.size(); i++){
		map<int, int>::iterator it = ans.find((*(data[i]))[targetCol]);
		if(it != ans.end()){
			ans[(*(data[i]))[targetCol]] += 1;
		}
		else{
			ans[(*(data[i]))[targetCol]] = 1;
		}
	}
	int maxNumb=0, maxClass = 0;
	for(map<int, int>::iterator it = ans.begin(); it!=ans.end(); it++){
		if(maxNumb < it->second){
			maxNumb = it->second;
			maxClass = it->first;
		}
	}
	return maxClass;
}

//检查节点所含类个数 大于1返回false 为1返回true
bool check_class_numb(vector<vector<int>* > data){

    set<int> classNumb;
    for(int iData=0; iData<data.size(); iData++){
        classNumb.insert((*(data[iData]))[targetCol]);
    }
    return classNumb.size()>1 ? false : true;
}

//计算GINI指数
double get_gini_val(vector<vector<int>* > data, int featureTemp){
    double gini = 0.0;

    for(int FiLoop=0; FiLoop<d_feature_numb[featureTemp].size(); FiLoop++){

        //计算特征值为Fi 有多少数据
        int Fi = d_feature_numb[featureTemp][FiLoop];
        int dataFiNumb = 0;
        for(int iData=0; iData<data.size(); iData++)
            if((*(data[iData]))[featureTemp] == Fi)
                dataFiNumb++;
        if(dataFiNumb == 0) continue;

        //找到该子节点每个类所含数据数量
        map<int, int> classDataNumb;
        for(int iData=0; iData<data.size(); iData++){

            if((*(data[iData]))[featureTemp] == Fi) {
                map<int, int>::iterator classDataNumbIterator = classDataNumb.find((*(data[iData]))[targetCol]);
                if(classDataNumbIterator != classDataNumb.end())
                    classDataNumb[(*(data[iData]))[targetCol]] += 1;
                else
                    classDataNumb[(*(data[iData]))[targetCol]] = 1;
            }
        }

        double p = 0.0;
        for(map<int, int>::iterator classDataNumbIterator = classDataNumb.begin(); classDataNumbIterator!=classDataNumb.end(); classDataNumbIterator++){
            p += pow(classDataNumbIterator->second*1.0 / dataFiNumb, 2);
        }
        gini += dataFiNumb*1.0/data.size()*(1 - p);
    }
    return gini;
}

//统计不一致数据量
int get_diff_data_num(vector<vector<int>* > data){

    int diffNumb = 0;
    for(map<int, int>::iterator AIterator=d_kv_fd.begin(); AIterator!=d_kv_fd.end(); AIterator++){
        for(int AiLoop=0; AiLoop<d_feature_numb[AIterator->first].size(); AiLoop++){

            int AiDiffNumb = 0;
            set<int> BiDiffNumb;
            for(int iData=0; iData<data.size(); iData++){
                if((*(data[iData]))[AIterator->first] == d_feature_numb[AIterator->first][AiLoop]){
                    AiDiffNumb++;
                    BiDiffNumb.insert((*(data[iData]))[AIterator->second]);
                }
            }
            if(BiDiffNumb.size() > 1)
                diffNumb += AiDiffNumb;
        }
    }
    return diffNumb;
}

//计算分歧率
double get_fenqilv(vector<vector<int>* > data, int featureTemp, int diffNumb){
    double fenqilv = 0.0;
    for(int AiLoop=0; AiLoop<d_feature_numb[featureTemp].size(); AiLoop++){

        double Da = 0.0, fenqilvTemp = 1.0;
        map<int, int> BiNumb;
        for(int iData=0; iData<data.size(); iData++){
            if((*(data[iData]))[featureTemp] == d_feature_numb[featureTemp][AiLoop]){
                Da++;
                map<int, int>::iterator BiNumbIterator = BiNumb.find((*(data[iData]))[d_kv_fd[featureTemp]]);
                if(BiNumbIterator != BiNumb.end())
                    BiNumb[(*(data[iData]))[d_kv_fd[featureTemp]]]++;
                else
                    BiNumb[(*(data[iData]))[d_kv_fd[featureTemp]]] = 1;
            }
        }
        fenqilvTemp = Da*1.0/diffNumb;
        for(map<int, int>::iterator BiNumbIterator = BiNumb.begin(); BiNumbIterator!=BiNumb.end(); BiNumbIterator++){
            fenqilvTemp *= (BiNumbIterator->second*1.0/Da);
        }

        fenqilv += (Da*1.0/data.size()*fenqilvTemp);
    }
    return 1.0 - fenqilv;
}

//计算占比率
double get_zhanbilv(vector<vector<int>* > data, int featureTemp){

    double zhanbilv = 0.0;
    for(int BiLoop=0; BiLoop<d_feature_numb[featureTemp].size(); BiLoop++){

        double zhanbilvTemp = 1.0, subNodeNumb = 0.0;

        for(int iData=0; iData<data.size(); iData++){
            if((*(data[iData]))[featureTemp] == d_feature_numb[featureTemp][BiLoop])
                subNodeNumb++;
        }

        for(int AiLoop=0; AiLoop<d_feature_numb[d_kv_fd_2[featureTemp]].size(); AiLoop++){

            double subNodeAiNumb = 0.0, fatherNodeAiNumb = 0.0;
            for(int iData=0; iData<data.size(); iData++){

                //计算子节点Ai数
                if((*(data[iData]))[featureTemp] == d_feature_numb[featureTemp][BiLoop] &&
                        ((*(data[iData]))[d_kv_fd_2[featureTemp]] == d_feature_numb[d_kv_fd_2[featureTemp]][AiLoop])){
                    subNodeAiNumb++;
                }

                //计算父节点Ai数
                if((*(data[iData]))[d_kv_fd_2[featureTemp]] == d_feature_numb[d_kv_fd_2[featureTemp]][AiLoop]){
                    fatherNodeAiNumb++;
                }
            }
            if(subNodeAiNumb != 0 && fatherNodeAiNumb != 0)
                zhanbilvTemp *= (subNodeAiNumb*1.0/fatherNodeAiNumb);
        }
        zhanbilv += (subNodeNumb*1.0/data.size()*zhanbilvTemp);
    }
    return 1.0 - zhanbilv;
}

//计算占比率
double get_hunxiaolv(vector<vector<int>* > data, int featureTemp){
    double hunxiaolv = 0.0;
    for(int BiLoop=0; BiLoop<d_feature_numb[featureTemp].size(); BiLoop++){

        double hunxiaolvTemp = 1.0, subNodeNumb = 0.0;

        for(int iData=0; iData<data.size(); iData++){
            if((*(data[iData]))[featureTemp] == d_feature_numb[featureTemp][BiLoop])
                subNodeNumb++;
        }

        for(int AiLoop=0; AiLoop<d_feature_numb[d_kv_fd_2[featureTemp]].size(); AiLoop++){

            double subNodeAiNumb = 0.0;
            for(int iData=0; iData<data.size(); iData++){

                //计算子节点Ai数
                if((*(data[iData]))[featureTemp] == d_feature_numb[featureTemp][BiLoop] &&
                        ((*(data[iData]))[d_kv_fd_2[featureTemp]] == d_feature_numb[d_kv_fd_2[featureTemp]][AiLoop])){
                    subNodeAiNumb++;
                }
            }
            if(subNodeAiNumb != 0 && subNodeNumb != 0)
                hunxiaolvTemp *= (subNodeAiNumb*1.0/subNodeNumb);
        }
        hunxiaolv += (subNodeNumb*1.0/data.size()*hunxiaolvTemp);
    }
    return hunxiaolv!=1.0 ? hunxiaolv : 0.0;
}

Node* build_tree(int treeLevel, bool if_add_stra, bool if_print=true){
	Node* root = &(treeNode[treeNodeIndex++]);
	for(int iData=0; iData<train_data.size(); iData++){
		root->nodeData.push_back(&(train_data[iData]));
	}
	root->maxClass = set_class(root->nodeData);
    root->unique = treeNodeIndex-1;

	queue<Node*> nodeQueue;
	nodeQueue.push(root);
	while(!nodeQueue.empty()){

		map<int, double> giniList;
		Node* nodeTemp = nodeQueue.front();
		nodeQueue.pop();

        //如果该节点类个数为1 或者 该节点树层数大于指定层数 则不进行扩展
        if(check_class_numb(nodeTemp->nodeData) || nodeTemp->treeLevel >= treeLevel)
            continue;

		for(int iFeature=0; iFeature<l_feature.size(); iFeature++){

            double giniNow = 0.0;

            //分割的feature非函数依赖中的 或 是函数依赖中的前置但是对应的后置值唯一，只计算GINI指数
            map<int, int>::iterator kv1Iterator = d_kv_fd.find(l_feature[iFeature]);
            map<int, int>::iterator kv2Iterator = d_kv_fd_2.find(l_feature[iFeature]);

            if((kv1Iterator==d_kv_fd.end() && kv2Iterator==d_kv_fd_2.end()) ||
                    (kv1Iterator!=d_kv_fd.end() && d_feature_numb[l_feature[iFeature]].size()==1)){

                giniNow = get_gini_val(nodeTemp->nodeData, l_feature[iFeature]);

                if(if_print)
                    cout<<"GINI : "<<giniNow<<endl;
            }

            // 否则除计算GINI指数外，还需计算其他比率
            else{

                //feature是前置且对应的B未被分割过,需计算分歧率+GINI, 目标函数GINI+1-分歧率(Dai/总共不一致数据量 * Dbi/Dai)
                map<int, int>::iterator kvFiIterator = d_kv_fd.find(l_feature[iFeature]);
                set<int>::iterator usedIterator = l_used.find(d_kv_fd[l_feature[iFeature]]);
                map<int, int>::iterator kvFiIterator2 = d_kv_fd_2.find(l_feature[iFeature]);
                set<int>::iterator usedIterator2 = l_used.find(d_kv_fd_2[l_feature[iFeature]]);
                if(kvFiIterator != d_kv_fd.end() && usedIterator==l_used.end()){

                    double fenqilv = 0.0;
                    int diffDataNum = get_diff_data_num(nodeTemp->nodeData);

                    if(diffDataNum > 0){
                        fenqilv = get_fenqilv(nodeTemp->nodeData, l_feature[iFeature], diffDataNum);
                    }
                    giniNow = get_gini_val(nodeTemp->nodeData, l_feature[iFeature]);

                    if(if_print)
                        cout<<"GINI: "<<giniNow<<" fenqilv: "<<fenqilv<<endl;
                    if(if_add_stra)
                        giniNow += fenqilv;
                }
                //feature是后置且对应的A未被分割过,计算占比率+混淆率+GINI
                else if(kvFiIterator2 != d_kv_fd_2.end() && usedIterator2==l_used.end()){

                    double zhanbilv = 0.0, hunxiaolv = 0.0;
                    zhanbilv = get_zhanbilv(nodeTemp->nodeData, l_feature[iFeature]);
                    hunxiaolv = get_hunxiaolv(nodeTemp->nodeData, l_feature[iFeature]);
                    giniNow = get_gini_val(nodeTemp->nodeData, l_feature[iFeature]);

                    if(if_print)
                        cout<<"GINI: "<<giniNow<<" zhanbilv: "<<zhanbilv<<" huxiaolv: "<<hunxiaolv<<endl;
                    if(if_add_stra)
                        giniNow += (zhanbilv + hunxiaolv);
                } else
                    giniNow = get_gini_val(nodeTemp->nodeData, l_feature[iFeature]);
            }
            giniList[l_feature[iFeature]] = giniNow;

            if(if_print)
                cout<<"Ending of calculate feature: "<<l_feature[iFeature]<<" Gini: "<<giniNow<<" treeLevel "<<nodeTemp->treeLevel<<endl;
		}

        if(giniList.empty()) continue;
        int minFeature = 0;
        double minGini = (numeric_limits<double>::max)();
        for(map<int, double>::iterator giniIterator=giniList.begin(); giniIterator!=giniList.end(); giniIterator++){
            if(giniIterator->second < minGini){
                minGini = giniIterator->second;
                minFeature = giniIterator->first;
            }
        }
        if(if_print)
            cout<<"Best feature is: "<<minFeature<<" in all of feature, GINI' is: "<<minGini<<endl;
        l_used.insert(minFeature);

        nodeTemp->splitNum = minFeature;
        for(int FiLoop=0; FiLoop<d_feature_numb[minFeature].size(); FiLoop++){
            vector<vector<int>* > subFiNodeData;
            for(int iData=0; iData<nodeTemp->nodeData.size(); iData++){
                if((*(nodeTemp->nodeData[iData]))[minFeature] == d_feature_numb[minFeature][FiLoop])
                    subFiNodeData.push_back(nodeTemp->nodeData[iData]);
            }
            if(!subFiNodeData.empty()){
                Node* subNodeTemp = &(treeNode[treeNodeIndex++]);
                subNodeTemp->nodeData = subFiNodeData;
                subNodeTemp->maxClass = set_class(subFiNodeData);
                subNodeTemp->treeLevel = nodeTemp->treeLevel+1;
                subNodeTemp->unique = treeNodeIndex-1;
                (nodeTemp->subNode)[d_feature_numb[minFeature][FiLoop]] = subNodeTemp;
                nodeQueue.push(subNodeTemp);
            }
        }
	}
    return root;
}

// 深度复制参数树模型结构
Node* copy_tree(Node* root){
    Node* rootTemp = new Node(root->splitNum, root->maxClass, root->treeLevel, root->unique);
    rootTemp->subNode = root->subNode;
    rootTemp->nodeData = root->nodeData;

    queue<Node*> queueNode;
    queueNode.push(rootTemp);

    while(!queueNode.empty()){

        Node* temp = queueNode.front();
        queueNode.pop();

        map<int, Node*> subNodeNew;
        for(map<int, Node*>::iterator it=temp->subNode.begin(); it!=temp->subNode.end(); it++){

            Node* nodeTemp = new Node((*(it->second)).splitNum, (*(it->second)).maxClass, (*(it->second)).treeLevel, (*(it->second)).unique);
            nodeTemp->nodeData = (*(it->second)).nodeData;
            nodeTemp->subNode = (*(it->second)).subNode;

            subNodeNew[it->first] = nodeTemp;
            queueNode.push(nodeTemp);

        }
        temp->subNode = subNodeNew;
    }

    return rootTemp;
}

// 计算节点的不剪枝误差代价， 以及节点个数
void cul_node_err_sum_and_node_numb(Node root, double* rtt, double* nodeNumb){

    queue<Node> queueNode;
    queueNode.push(root);

    while(!queueNode.empty()){
        Node tempNode = queueNode.front();
        queueNode.pop();

        if(tempNode.subNode.empty()){
            int otherClassDataNumb = 0;
            for(vector<vector<int>* >::iterator dataIterator=tempNode.nodeData.begin(); dataIterator!=tempNode.nodeData.end(); dataIterator++)
                if((*(*dataIterator))[targetCol] != tempNode.maxClass)
                    otherClassDataNumb++;
            *rtt = (*rtt) + (otherClassDataNumb*1.0/root.nodeData.size());
            *nodeNumb = (*nodeNumb)+1;
        }
        else{
            for(map<int, Node*>::iterator subNodeIterator=tempNode.subNode.begin(); subNodeIterator!=tempNode.subNode.end(); subNodeIterator++)
                queueNode.push(*(subNodeIterator->second));
        }
    }
}

// 计算剪枝策略中的alpha值
double cul_pruning_alpha(Node root){

    int otherClassDataNumb = 0;
    for(vector<vector<int>* >::iterator dataIterator=root.nodeData.begin(); dataIterator!=root.nodeData.end(); dataIterator++)
        if ((*(*dataIterator))[targetCol] != root.maxClass)
            otherClassDataNumb++;

    double rt = otherClassDataNumb*1.0 / root.nodeData.size();
    double rtt = 0, node_numb = 0;
    cul_node_err_sum_and_node_numb(root, &rtt, &node_numb);
    double alpha = (rt - rtt) / (node_numb - 1.0);

    return alpha;
}

//剪枝, 每次只剪一个节点，返回经过alpha剪枝后的子树Node
Node* pruning_one(Node* root){

    Node* copyRoot = copy_tree(root);

    queue<Node*> queueNode;
    queueNode.push(copyRoot);
    //记录address(Node):alpha值  对每个节点剪枝一次 得到n+1个子树 从中取得最好的子树
    map<long long, double> alpha;

    //计算每个非叶子节点的alpha值
    while (!queueNode.empty()){

        Node temp = *(queueNode.front());
        queueNode.pop();

        if(!temp.subNode.empty())
            alpha[temp.unique] = cul_pruning_alpha(temp);
        for(map<int, Node*>::iterator subNodeIterator=temp.subNode.begin(); subNodeIterator!=temp.subNode.end(); subNodeIterator++)
            queueNode.push(subNodeIterator->second);
    }

    long long minNodeAdr = -2;
    double minNodeAlpha = (numeric_limits<double>::max)();
    for(map<long long, double>::iterator alphaIterator=alpha.begin(); alphaIterator!=alpha.end(); alphaIterator++){
        if(alphaIterator->second < minNodeAlpha){
            minNodeAlpha = alphaIterator->second;
            minNodeAdr = alphaIterator->first;
        }
    }

    //剪去最小alpha对应的非叶子节点
    queueNode.push(copyRoot);
    while(!queueNode.empty()){
        Node* temp = queueNode.front();
        queueNode.pop();
        if(!temp->subNode.empty()) {
            if (temp->unique == minNodeAdr) {
                temp->subNode.clear();
            }
            else {
                for (map<int, Node*>::iterator subNodeIterator = temp->subNode.begin();
                     subNodeIterator != temp->subNode.end(); subNodeIterator++)
                    queueNode.push(subNodeIterator->second);
            }
        }
    }
//    cout<<"before pruing node numb: "<<check_node_numb(*root)<<" after: "<<check_node_numb(*copyRoot)<<endl;
    return copyRoot;
}

//剪枝, 调用pruning_one每次剪出最优子树,不断剪枝直至只剩根节点，返回经过alpha剪枝后的子树列表集合
queue<Node*> pruning_all(Node* root){
    Node* nodeTemp = root;
    queue<Node*> queueTree;
    queueTree.push(nodeTemp);
    while (!nodeTemp->subNode.empty()){
        Node* subTree = pruning_one(nodeTemp);
        queueTree.push(subTree);
        nodeTemp = subTree;
    }
    return queueTree;
}

// 分类
double classify(Node* root){

    int errDataNumb = 0;

    for(int iData=0; iData<test_data.size(); iData++){
        Node* temp = root;
        while(temp->splitNum != -1){
            bool findAns = false;
            for(map<int, Node*>::iterator subNodeIterator=temp->subNode.begin(); subNodeIterator!=temp->subNode.end(); subNodeIterator++){
                if(test_data[iData][temp->splitNum] == subNodeIterator->first) {
                    temp = temp->subNode[subNodeIterator->first];
                    findAns = true;
                    break;
                }
            }
            if(!findAns)
                break;
        }
        if(temp->maxClass != test_data[iData][targetCol])
            errDataNumb++;
    }

    return errDataNumb*1.0 / test_data.size();
}

//void delete_tree(Node* root){
//    if(root->subNode.empty()){
//        delete root;
//        return;
//    }
//    for(map<int, Node*>::iterator subNodeIterator=root->subNode.begin(); subNodeIterator!=root->subNode.end(); subNodeIterator++)
//        delete_tree(subNodeIterator->second);
//    return;
//}

// 从子树列表中，根据测试集，选出最优的子树
void choose_best_tree(Node* root, queue<Node*> sub_tree_queue){

    double minErr = 100.0;

    while(!sub_tree_queue.empty()){
        Node* tempTree = sub_tree_queue.front();
        sub_tree_queue.pop();
        double errTemp = classify(tempTree);
        if(errTemp < minErr){
            minErr = errTemp;
        }
//        if(root != tempTree)
//            delete_tree(tempTree);
    }

    cout<<"ANS: "<<1.0-minErr<<endl;
    return ;
}

void pruning(Node* root){
    queue<Node*> sub_tree_queue = pruning_all(root);
//    cout<<"sub_tree_queue.size "<<sub_tree_queue.size()<<endl;
    choose_best_tree(root, sub_tree_queue);
    return ;
}


int main(){

    train_data.clear();
    test_data.clear();
    treeNodeIndex = 0;
    data_max_size = 5000;
    int treeLevel = 100;
    int i = 6;
    string adr = "C:\\Users\\Whiker\\Desktop\\Tree_C++\\krkopt-changed-";
    char buffer[10];
    adr = adr + string(_itoa(i, buffer, 10)) + ".csv";
//    string adr = "C:\\Users\\Whiker\\Desktop\\Tree_C++\\krkopt-changed.csv";

    vector<vector<int> > data = read_data(adr);
    pre_data(data);
    cout << i << " per data size: " << data.size() << " train data size: " << train_data.size()
         << " test data size: " << test_data.size() << endl;

    Node *root = build_tree(treeLevel, true, false);
//    cout<<"build tree root numb: "<<check_node_numb(root)<<endl;

//    cout<<i<<" per  add stra "<<1.0-classify(root)<<endl;
    pruning(root);
//        Node bestTree = pruning(root);
//        cout<<i<<" per  !!!!!!!!!!!!!!!add stra ANS: "<<1.0-classify(&bestTree)<<endl;

    Node *rootNoStra = build_tree(treeLevel, false, false);
//    cout<<i<<" per  no add stra "<<1.0-classify(rootNoStra)<<endl;
    pruning(rootNoStra);
//        Node bestTreeNoStra = pruning(rootNoStra);
//        cout<<i<<" per  !!!!!!!!!!!!!!!no stra ANS: "<<1.0-classify(&bestTreeNoStra)<<endl;
    return 0;
}

/*
//验证 d_kv_fd_numb  d_feature_numb  d_maxBi准确性
void yanzheng_predata(){
	cout<<"d_kv_fd_numb :"<<endl;
	for(map<int, map<int, map<int,int> > >::iterator i=d_kv_fd_numb.begin(); i!=d_kv_fd_numb.end(); i++){
		cout<<"A: "<<i->first<<endl;
		for(map<int, map<int,int> >::iterator ii=i->second.begin(); ii!=i->second.end(); ii++){
			cout<<"Ai: "<<ii->first<<endl;
			for(map<int, int>::iterator iii=ii->second.begin(); iii!=ii->second.end(); iii++){
				cout<<"Bi:numb "<<iii->first<<":"<<iii->second<<endl;
			}
		}
		cout<<endl;
	}

	cout<<"d_feature_numb :"<<endl;
	for(map<int, vector<int> >::iterator i=d_feature_numb.begin(); i!=d_feature_numb.end(); i++){
		cout<<endl<<"A: "<<i->first<<endl<<"Ai:";
		for(int ii=0; ii<i->second.size(); ii++){
			cout<<"   "<<i->second[ii];
		}
		cout<<endl;
	}

	cout<<"d_maxBi :"<<endl;
	for(map<int, map<int, int> >::iterator i=d_maxBi.begin(); i!=d_maxBi.end(); i++){
		cout<<"A: "<<i->first<<endl;
		for(map<int, int>::iterator ii=i->second.begin(); ii!=i->second.end(); ii++){
			cout<<"Ai:Bi "<<ii->first<<":"<<ii->second<<endl;
		}
		cout<<endl;
	}
}

long long yanzhengIndex = 0;
void yanzheng_tree(Node root){
    yanzhengIndex++;
    if(!root.subNode.empty()){
        for(int i=0; i<root.subNode.size(); i++)
            yanzheng_tree(*(root.subNode[i]));
    }
}
//检查两个树模型地址值是否一致
void check_copy(Node* root, Node* copyRoot){
    queue<Node*> queue1;
    queue<Node*> queue2;
    queue1.push(root);
    queue2.push(copyRoot);
    long long sameNumb = 0;
    while(!queue1.empty()){
        Node* temp1 = queue1.front();
        queue1.pop();
        Node* temp2 = queue2.front();
        queue2.pop();
        if(temp1 == temp2){
            sameNumb++;
        }
        if(!(temp1->subNode.empty())){
            queue1.push(temp1->subNode[temp1->splitFi]);
            queue1.push(temp1->subNode[temp1->splitFi + 1]);
        }
        if(!(temp2->subNode.empty())){
            queue2.push(temp2->subNode[temp2->splitFi]);
            queue2.push(temp2->subNode[temp2->splitFi + 1]);
        }
    }
    cout<<"same numb: "<<sameNumb<<endl;
}
*/