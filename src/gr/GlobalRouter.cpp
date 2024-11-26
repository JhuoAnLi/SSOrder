#include "GlobalRouter.h"
#include "PatternRoute.h"
#include "MazeRoute.h"
#include <chrono>
#include <queue>
#include <thread>
#include <mutex>
#include <future>


GlobalRouter::GlobalRouter(const Design& design, const Parameters& params): 
    gridGraph(design, params), parameters(params) {
    // Instantiate the global routing netlist
    const vector<Net>& baseNets = design.getAllNets();
    nets.reserve(baseNets.size());
    for (const Net& baseNet : baseNets) {
        nets.emplace_back(baseNet, design, gridGraph);
    }
}

std :: priority_queue<std :: pair<double,unsigned>> GlobalRouter::selectNext(const PatternRoute& patternRoute, GridGraphView<CostT>& rudyView, const double& AverageResource, const vector<unsigned>& ShortNewnetIndices, vector<vector<CapacityT>>& TdCapacity) {
        std::priority_queue<std::pair<double, unsigned>> pq;
        std::vector<std::vector<int>> sameDirectionLayers(2);
        double TotalBeRoutedNet(0.0);

        for (const auto& row : TdCapacity) {
            for (const auto& capacity : row) {
                TotalBeRoutedNet += (capacity - AverageResource) * (capacity - AverageResource);
            }
        }
        for (unsigned layerIndex = parameters.min_routing_layer; layerIndex < gridGraph.getNumLayers(); layerIndex++)
            sameDirectionLayers[gridGraph.getLayerDirection(layerIndex)].emplace_back(layerIndex);
        std::vector<std::pair<utils::PointT<int>, utils::PointT<int>>> edges;
        SteinerTreeNode::preorder(patternRoute.steinerTree, [&edges](std::shared_ptr<SteinerTreeNode> node) {
            for (auto& child : node->children) {
                edges.emplace_back(*node, *child);
            }
        });

        std::vector<std::future<void>> futures;
        std::mutex pq_mutex;
        unsigned n = ShortNewnetIndices.size();
        unsigned num_threads = std::min((unsigned)8, n);
        unsigned chunk_size = (ShortNewnetIndices.size() + num_threads - 1) / num_threads;

        double GridArea = gridGraph.getSize(0) * gridGraph.getSize(1);

        for (unsigned thread_id = 0; thread_id < num_threads; ++thread_id) {
            unsigned start_idx = thread_id * chunk_size;
            unsigned end_idx = std::min(start_idx + chunk_size, n);

            futures.emplace_back(std::async(std::launch::async, [&, start_idx, end_idx] {
                for (unsigned i = start_idx; i < end_idx; ++i) {
                    auto& net(nets[ShortNewnetIndices[i]]);
                    double NAR_Resource = gridGraph.lCost(sameDirectionLayers, edges, rudyView, net);

                    double LocalTotalBeRoutedNet = TotalBeRoutedNet;
                    double LedgeCost(0.0);
                    
                    auto bbox(net.getBoundingBox());
                    unsigned lx(bbox.lx()), hx(bbox.hx()), ly(bbox.ly()), hy(bbox.hy());

                    for (unsigned j = lx; j <= hx; ++j) {
                        // LedgeCost += 0.5 * (2 * TdCapacity[j][ly] - 2 * AverageResource - 0.5) + 0.5 * (2 * TdCapacity[j][hy] - 2 * AverageResource - 0.5);
                        // LedgeCost += 0.5 * ((2 * TdCapacity[j][ly] - 2 * AverageResource - 0.5) + (2 * TdCapacity[j][hy] - 2 * AverageResource - 0.5));
                        // LedgeCost += 0.5 * (2 * TdCapacity[j][ly] + 2 * TdCapacity[j][hy] - 4 * AverageResource - 1);
                        LedgeCost += TdCapacity[j][ly] + TdCapacity[j][hy] - 2 * AverageResource - 0.5;
                    }
                    for (unsigned j = ly+1; j < hy; ++j) {
                        LedgeCost += TdCapacity[lx][j] + TdCapacity[hx][j] - 2 * AverageResource - 0.5;
                    }

                    LocalTotalBeRoutedNet -= LedgeCost;
                    double VarianceSource = LocalTotalBeRoutedNet / GridArea;

                    std::lock_guard<std::mutex> lock(pq_mutex);
                    pq.push({NAR_Resource - VarianceSource, ShortNewnetIndices[i]});
                }
            }));
        }
        for (auto& future : futures) {
            future.get();
        }
        return pq;
    }

void GlobalRouter::route() {
    int n1 = 0, n2 = 0, n3 = 0;
    double t1 = 0, t2 = 0, t3 = 0;
    
    auto t = std::chrono::high_resolution_clock::now();
    
    vector<int> netIndices;
    netIndices.reserve(nets.size());
    for (const auto& net : nets) netIndices.push_back(net.getIndex());
    // Stage 1: Pattern routing
    n1 = netIndices.size();
    PatternRoute::readFluteLUT();
    log() << "stage 1: pattern routing" << std::endl;
    sortNetIndices(netIndices);

    GridGraphView<CostT> rudyView; // (2d) direction -> x -> y -> rudy value
    gridGraph.extractRUDY(rudyView, nets);
    GridGraphView<CostT> rudyView_bk = rudyView;
    // sortNetIndicesByRudy(netIndices, rudyView_bk);
    double acc_select = 0.0;
    double acc_construct = 0.0;
    double acc_run = 0.0;
    double acc_mse = 0.0;
    // log() << "Monotonic routing..." << std::endl;

    //#### JALI ####
    // vector<int> NewnetIndices = netIndices;
    // unsigned netIndex = NewnetIndices[0];
    // auto it = find(NewnetIndices.begin(), NewnetIndices.end(), netIndex);
    // double CurrentHp = nets[netIndex].getBoundingBox().hp();
    // unsigned netCnt = 0;
    // std :: priority_queue<std :: pair<double,unsigned>> pq;
    // vector<unsigned> ShortNewnetIndices;

    // double TotalCapacity = 0.0;
    // vector<vector<CapacityT>> TdCapacity(gridGraph.getSize(0), vector<CapacityT>(gridGraph.getSize(1), 0.0));
    //     for (int layerIndex = parameters.min_routing_layer; layerIndex < gridGraph.getNumLayers(); layerIndex++) {
    //     unsigned direction = gridGraph.getLayerDirection(layerIndex);
    //         for (int x = 0; x < gridGraph.getSize(0) - 1 + direction; x++) {
    //             for (int y = 0; y < gridGraph.getSize(1) - direction; y++) {
    //                 CapacityT capacity = max(gridGraph.getEdge(layerIndex, x, y).capacity, 0.0);
    //                 TotalCapacity += capacity;
    //                 TdCapacity[x][y] += capacity;
    //             }
    //         }
    //     }
    
    // do {
    //     PatternRoute patternRoute(nets[netIndex], gridGraph, parameters);
    //     patternRoute.constructSteinerTree();
    //     auto bbox(nets[netIndex].getBoundingBox());
    //     auto start = clock();
    //     // if((bbox.lx() != bbox.hx()) && (bbox.ly() != bbox.hy()))
    //     //     patternRoute.select(rudyView);
    //     auto end = clock();
    //     acc_select += (double(end-start)/CLOCKS_PER_SEC);
    //     start = clock();
    //     patternRoute.constructRoutingDAG();
    //     end = clock();
    //     acc_construct += (double(end-start)/CLOCKS_PER_SEC);
    //     patternRoute.run();
    //     start = clock();
    //     acc_run += (double(start - end)/CLOCKS_PER_SEC);
    //     gridGraph.commitTree(nets[netIndex].getRoutingTree());
    //     gridGraph.updateRUDY(rudyView, nets[netIndex]);

    //     if(pq.empty()) 
    //     {
    //         // gridGraph.write("./heatmap_temp/18t10m/heatmap_dowhile_"+ std::to_string(netCnt) +".txt");
    //         for (auto i = it + 1; i != NewnetIndices.end() && i < it + 1001 ; ++i) {
    //             ShortNewnetIndices.push_back(*i);
    //         }

    //         if(!ShortNewnetIndices.empty())
    //         {
    //             // log() << "SELECTING NEXT NET..."<< std::endl;

    //             it = find(NewnetIndices.begin(),NewnetIndices.end(),ShortNewnetIndices[ShortNewnetIndices.size()-1]);
    //             // CurrentHp += nets[netIndex].getBoundingBox().hp();
    //             double AverageResource((TotalCapacity - CurrentHp) / (gridGraph.getSize(0) * gridGraph.getSize(1)));
    //             start = clock();
    //             pq = selectNext(patternRoute, rudyView, AverageResource, ShortNewnetIndices, TdCapacity);
    //             end = clock();
    //             acc_mse += (double(end-start)/CLOCKS_PER_SEC);
    //             for(unsigned i = 0; i<ShortNewnetIndices.size();++i)
    //                 CurrentHp += nets[ShortNewnetIndices[i]].getBoundingBox().hp();
    //             ShortNewnetIndices.clear();
    //             // log() << "ROUTING..." << std::endl;
    //         }
    //     }
    //     if(!pq.empty())
    //     {
    //         netIndex = pq.top().second;
    //         pq.pop();
    //     }
    //     netCnt++;
    // } while (netCnt < netIndices.size());
    // gridGraph.write("./heatmap_jali/18t10m/heatmap_final.txt");
    // #### JALI ####
    for (const int netIndex : netIndices) {
        PatternRoute patternRoute(nets[netIndex], gridGraph, parameters);
        patternRoute.constructSteinerTree();
        
        auto bbox = nets[netIndex].getBoundingBox();
        auto start = clock();
        if((bbox.lx() != bbox.hx()) && (bbox.ly() != bbox.hy()))
            patternRoute.select(rudyView);
        auto end = clock();
        acc_select += (double(end-start)/CLOCKS_PER_SEC);

        start = clock();
        patternRoute.constructRoutingDAG();
        end = clock();
        acc_construct += (double(end-start)/CLOCKS_PER_SEC);
        patternRoute.run();
        start = clock();
        acc_run += (double(start - end)/CLOCKS_PER_SEC);
        gridGraph.commitTree(nets[netIndex].getRoutingTree());
        gridGraph.updateRUDY(rudyView, nets[netIndex]);
    }
    writeGain();
    
    std::cout << "[RUNTIME] select: " << acc_select << "s. construct: " << acc_construct << "s. run: " << acc_run << "s. mse_function: " << acc_mse <<"s."<<std::endl;
    // gridGraph.write("./heatmap_jali/18t10m_cugr2/heatmap_final.txt");
    netIndices.clear();
    for (const auto& net : nets) {
        if (gridGraph.checkOverflow(net.getRoutingTree()) > 0) {
            netIndices.push_back(net.getIndex());
        }
    }
    log() << netIndices.size() << " / " << nets.size() << " nets have overflows." << std::endl;
    logeol();
    
    t1 = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t).count();
    t = std::chrono::high_resolution_clock::now();
    
    // Stage 2: Pattern routing with possible detours
    
    n2 = netIndices.size();
    if (netIndices.size() > 0) {
        log() << "stage 2: pattern routing with possible detours" << std::endl;
        GridGraphView<bool> congestionView; // (2d) direction -> x -> y -> has overflow?
        gridGraph.extractCongestionView(congestionView);
        sortNetIndices(netIndices);
        // sortNetIndicesByCongest(netIndices, congestionView);
        for (const int netIndex : netIndices) {
            GRNet& net = nets[netIndex];
            net.setChoice(0);
            gridGraph.commitTree(net.getRoutingTree(), true);
            PatternRoute patternRoute(net, gridGraph, parameters);
            patternRoute.constructSteinerTree();
    
            patternRoute.constructRoutingDAG();
            if(net.getChoice() == 0) {
                patternRoute.constructDetours(congestionView); // KEY DIFFERENCE compared to stage 1
            }
            patternRoute.run();
            gridGraph.commitTree(net.getRoutingTree());
        }
        
        netIndices.clear();
        for (const auto& net : nets) {
            if (gridGraph.checkOverflow(net.getRoutingTree()) > 0) {
                netIndices.push_back(net.getIndex());
            }
        }
        log() << netIndices.size() << " / " << nets.size() << " nets have overflows." << std::endl;
        logeol();
    }
    // gridGraph.write("./heatmap_temp/18t10m/heatmap_2.txt");
    t2 = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t).count();
    t = std::chrono::high_resolution_clock::now();
    
    // Stage 3: maze routing on sparsified routing graph
    n3 = netIndices.size();
    if (netIndices.size() > 0) {
        log() << "stage 3: maze routing on sparsified routing graph" << std::endl;
        for (const int netIndex : netIndices) {
            GRNet& net = nets[netIndex];
            gridGraph.commitTree(net.getRoutingTree(), true);
        }
        GridGraphView<CostT> wireCostView;
        gridGraph.extractWireCostView(wireCostView);
        sortNetIndices(netIndices);
        // sortNetIndicesByRudy(netIndices, wireCostView);
        SparseGrid grid(10, 10, 0, 0);
        for (const int netIndex : netIndices) {
            GRNet& net = nets[netIndex];
            // gridGraph.commitTree(net.getRoutingTree(), true);
            // gridGraph.updateWireCostView(wireCostView, net.getRoutingTree());
            MazeRoute mazeRoute(net, gridGraph, parameters);
            mazeRoute.constructSparsifiedGraph(wireCostView, grid);
            mazeRoute.run();
            std::shared_ptr<SteinerTreeNode> tree = mazeRoute.getSteinerTree();
            assert(tree != nullptr);
            
            PatternRoute patternRoute(net, gridGraph, parameters);
            patternRoute.setSteinerTree(tree);
            patternRoute.constructRoutingDAG();
            patternRoute.run();
            
            gridGraph.commitTree(net.getRoutingTree());
            gridGraph.updateWireCostView(wireCostView, net.getRoutingTree());
            grid.step();
        }
        netIndices.clear();
        for (const auto& net : nets) {
            if (gridGraph.checkOverflow(net.getRoutingTree()) > 0) {
                netIndices.push_back(net.getIndex());
            }
        }
        log() << netIndices.size() << " / " << nets.size() << " nets have overflows." << std::endl;
        logeol();
    }
    // gridGraph.write("./heatmap_temp/18t10m/heatmap_3.txt");
    
    t3 = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t).count();
    t = std::chrono::high_resolution_clock::now();
    
    // std::cout << "iteration statistics " 
    //     << n1 << " " << std::setprecision(3) << std::fixed << t1 << " " 
    //     << n2 << " " << std::setprecision(3) << std::fixed << t2 << " " 
    //     << n3 << " " << std::setprecision(3) << std::fixed << t3 << std::endl;
    
    printStatistics();
    if (parameters.write_heatmap) gridGraph.write();
}

// void GlobalRouter::sortNetIndices(vector<int>& netIndices) const {
//     vector<int> halfParameters(nets.size());
    
//     for (int netIndex : netIndices) {
//         auto& net = nets[netIndex];
//         halfParameters[netIndex] = net.getBoundingBox().hp();
//     }

//     sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
//         return halfParameters[lhs] > halfParameters[rhs];
//     });
// }

void GlobalRouter::sortNetIndices(vector<int>& netIndices) const {
    vector<int> halfParameters(nets.size());
    
    for (int netIndex : netIndices) {
        auto& net = nets[netIndex];
        halfParameters[netIndex] = net.getBoundingBox().hp();
    }

    sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
        return halfParameters[lhs] < halfParameters[rhs];
    });
}

// void GlobalRouter::sortNetIndices(vector<int>& netIndices) const {
//     vector<double> halfParameters(nets.size());
//     //sort by area
//     for (int netIndex : netIndices) {
//         auto& net = nets[netIndex].getBoundingBox();
//         halfParameters[netIndex] = (net.hx() - net.lx()) * (net.hy() - net.ly());
//     }

//     sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
//         return halfParameters[lhs] > halfParameters[rhs];
//     });
// }

void GlobalRouter::sortNetIndicesByRudy(vector<int>& netIndices, GridGraphView<CostT>& rudy) const {
    vector<CostT> rudy_value(nets.size(), 0.0);
    
    for (int netIndex : netIndices) {
        auto& net = nets[netIndex];
        for(int x = net.getBoundingBox().lx(); x <= net.getBoundingBox().hx(); x++) {
            for(int y = net.getBoundingBox().ly(); y <= net.getBoundingBox().hy(); y++) {
                rudy_value[netIndex] += rudy[MetalLayer::H][x][y];
                rudy_value[netIndex] += rudy[MetalLayer::V][x][y];
            }
        }
        
        rudy_value[netIndex] /= net.getBoundingBox().hp();
    }

    sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
        return rudy_value[lhs] > rudy_value[rhs];
    });
}

void GlobalRouter::sortNetIndicesByCongest(vector<int>& netIndices, GridGraphView<bool>& congest) const {
    vector<CostT> congest_value(nets.size(), 0.0);
    
    for (int netIndex : netIndices) {
        auto& net = nets[netIndex];
        for(int x = net.getBoundingBox().lx(); x <= net.getBoundingBox().hx(); x++) {
            for(int y = net.getBoundingBox().ly(); y <= net.getBoundingBox().hy(); y++) {
                congest_value[netIndex] += double(congest[MetalLayer::H][x][y]);
                congest_value[netIndex] += double(congest[MetalLayer::V][x][y]);
            }
        }
        congest_value[netIndex] /= net.getBoundingBox().hp();
    }

    sort(netIndices.begin(), netIndices.end(), [&](int lhs, int rhs) {
        return congest_value[lhs] > congest_value[rhs];
    });
}

void GlobalRouter::getGuides(const GRNet& net, vector<std::pair<int, utils::BoxT<int>>>& guides) {
    auto& routingTree = net.getRoutingTree();
    if (!routingTree) return;
    // 0. Basic guides
    GRTreeNode::preorder(routingTree, [&](std::shared_ptr<GRTreeNode> node) {
        for (const auto& child : node->children) {
            if (node->layerIdx == child->layerIdx) {
                guides.emplace_back(
                    node->layerIdx, utils::BoxT<int>(
                        min(node->x, child->x), min(node->y, child->y),
                        max(node->x, child->x), max(node->y, child->y)
                    )
                );
            } else {
                int maxLayerIndex = max(node->layerIdx, child->layerIdx);
                for (int layerIdx = min(node->layerIdx, child->layerIdx); layerIdx <= maxLayerIndex; layerIdx++) {
                    guides.emplace_back(layerIdx, utils::BoxT<int>(node->x, node->y));
                }
            }
        }
    });
    
    
    auto getSpareResource = [&] (const GRPoint& point) {
        double resource = std::numeric_limits<double>::max();
        unsigned direction = gridGraph.getLayerDirection(point.layerIdx);
        if (point[direction] + 1 < gridGraph.getSize(direction)) {
            resource = min(resource, gridGraph.getEdge(point.layerIdx, point.x, point.y).getResource());
        }
        if (point[direction] > 0) {
            GRPoint lower = point;
            lower[direction] -= 1;
            resource = min(resource, gridGraph.getEdge(lower.layerIdx, point.x, point.y).getResource());
        }
        return resource;
    };
    
    // 1. Pin access patches
    assert(parameters.min_routing_layer + 1 < gridGraph.getNumLayers());
    for (auto& gpts : net.getPinAccessPoints()) {
        for (auto& gpt : gpts) {
            if (gpt.layerIdx < parameters.min_routing_layer) {
                int padding = 0;
                if (getSpareResource({parameters.min_routing_layer, gpt.x, gpt.y}) < parameters.pin_patch_threshold) {
                    padding = parameters.pin_patch_padding;
                }
                for (int layerIdx = gpt.layerIdx; layerIdx <= parameters.min_routing_layer + 1; layerIdx++) {
                    guides.emplace_back(layerIdx, utils::BoxT<int>(
                        max(gpt.x - padding, 0),
                        max(gpt.y - padding, 0),
                        min(gpt.x + padding, (int)gridGraph.getSize(0) - 1),
                        min(gpt.y + padding, (int)gridGraph.getSize(1) - 1)
                    ));
                    areaOfPinPatches += (guides.back().second.x.range() + 1) * (guides.back().second.y.range() + 1);
                }
            }
        }
    }
    
    // 2. Wire segment patches
    GRTreeNode::preorder(routingTree, [&](std::shared_ptr<GRTreeNode> node) {
        for (const auto& child : node->children) {
            if (node->layerIdx == child->layerIdx) {
                double wire_patch_threshold = parameters.wire_patch_threshold;
                unsigned direction = gridGraph.getLayerDirection(node->layerIdx);
                int l = min((*node)[direction], (*child)[direction]);
                int h = max((*node)[direction], (*child)[direction]);
                int r = (*node)[1 - direction];
                for (int c = l; c <= h; c++) {
                    bool patched = false;
                    GRPoint point = (direction == MetalLayer::H ? GRPoint(node->layerIdx, c, r) : GRPoint(node->layerIdx, r, c));
                    if (getSpareResource(point) < wire_patch_threshold) {
                        for (int layerIndex = node->layerIdx - 1; layerIndex <= node->layerIdx + 1; layerIndex += 2) {
                            if (layerIndex < parameters.min_routing_layer || layerIndex >= gridGraph.getNumLayers()) continue;
                            if (getSpareResource({layerIndex, point.x, point.y}) >= 1.0) {
                                guides.emplace_back(layerIndex, utils::BoxT<int>(point.x, point.y));
                                areaOfWirePatches += 1;
                                patched = true;
                            }
                        }
                    } 
                    if (patched) {
                        wire_patch_threshold = parameters.wire_patch_threshold;
                    } else {
                        wire_patch_threshold *= parameters.wire_patch_inflation_rate;
                    }
                }
            }
        }
    });
}

void GlobalRouter::printStatistics() const {
    log() << "routing statistics" << std::endl;
    loghline();

    // wire length and via count
    uint64_t wireLength = 0;
    int viaCount = 0;
    vector<vector<vector<int>>> wireUsage;
    wireUsage.assign(
        gridGraph.getNumLayers(), vector<vector<int>>(gridGraph.getSize(0), vector<int>(gridGraph.getSize(1), 0))
    );
    for (const auto& net : nets) {
        GRTreeNode::preorder(net.getRoutingTree(), [&] (std::shared_ptr<GRTreeNode> node) {
            for (const auto& child : node->children) {
                if (node->layerIdx == child->layerIdx) {
                    unsigned direction = gridGraph.getLayerDirection(node->layerIdx);
                    int l = min((*node)[direction], (*child)[direction]);
                    int h = max((*node)[direction], (*child)[direction]);
                    int r = (*node)[1 - direction];
                    for (int c = l; c < h; c++) {
                        wireLength += gridGraph.getEdgeLength(direction, c);
                        int x = direction == MetalLayer::H ? c : r;
                        int y = direction == MetalLayer::H ? r : c;
                        wireUsage[node->layerIdx][x][y] += 1;
                    }
                } else {
                    viaCount += abs(node->layerIdx - child->layerIdx);
                }
            }
        });
    }
    
    // resource
    CapacityT overflow = 0;

    CapacityT minResource = std::numeric_limits<CapacityT>::max();
    GRPoint bottleneck(-1, -1, -1);
    for (int layerIndex = parameters.min_routing_layer; layerIndex < gridGraph.getNumLayers(); layerIndex++) {
        unsigned direction = gridGraph.getLayerDirection(layerIndex);
        for (int x = 0; x < gridGraph.getSize(0) - 1 + direction; x++) {
            for (int y = 0; y < gridGraph.getSize(1) - direction; y++) {
                CapacityT resource = gridGraph.getEdge(layerIndex, x, y).getResource();
                if (resource < minResource) {
                    minResource = resource;
                    bottleneck = {layerIndex, x, y};
                }
                CapacityT usage = wireUsage[layerIndex][x][y];
                CapacityT capacity = max(gridGraph.getEdge(layerIndex, x, y).capacity, 0.0);
                if (usage > 0.0 && usage > capacity) {
                    overflow += usage - capacity;
                }
            }
        }
    }
    
    log() << "wire length (metric):  " << wireLength / gridGraph.getM2Pitch() << std::endl;
    log() << "total via count:       " << viaCount << std::endl;
    log() << "total wire overflow:   " << (int)overflow << std::endl;
    logeol();

    log() << "min resource: " << minResource << std::endl;
    log() << "bottleneck:   " << bottleneck << std::endl;

    logeol();
}

void GlobalRouter::writeGain(const std::string gain_file) const {
    log() << "writing gain to file..." << std::endl;
    std::stringstream ss;

    assert(vec_lGain.size() == vec_zGain.size() && vec_lGain.size() == vec_monoGain.size());
    for(auto lgain : vec_lGain){
        ss << lgain << " ";
    }
    ss << std::endl;

    for(auto zgain : vec_zGain){
        ss << zgain << " ";
    }
    ss << std::endl;

    for(auto monogain : vec_monoGain){
        ss << monogain << " ";
    }
    std::ofstream fout(gain_file);
    fout << ss.str();
    fout.close();
}

void GlobalRouter::write(std::string guide_file) {
    log() << "generating route guides..." << std::endl;
    if (guide_file == "") guide_file = parameters.out_file;
    
    areaOfPinPatches = 0;
    areaOfWirePatches = 0;
    std::stringstream ss;
    for (const GRNet& net : nets) {
        vector<std::pair<int, utils::BoxT<int>>> guides;
        getGuides(net, guides);
        
        ss << net.getName() << std::endl;
        ss << "(" << std::endl;
        for (const auto& guide : guides) {
            ss << gridGraph.getGridline(0, guide.second.x.low) << " "
                 << gridGraph.getGridline(1, guide.second.y.low) << " "
                 << gridGraph.getGridline(0, guide.second.x.high + 1) << " "
                 << gridGraph.getGridline(1, guide.second.y.high + 1) << " "
                 << gridGraph.getLayerName(guide.first) << std::endl;
        }
        ss << ")" << std::endl;
    }
    log() << "total area of pin access patches: " << areaOfPinPatches << std::endl;
    log() << "total area of wire segment patches: " << areaOfWirePatches << std::endl;
    log() << std::endl;
    log() << "writing output..." << std::endl;
    std::ofstream fout(guide_file);
    fout << ss.str();
    fout.close();
    log() << "finished writing output..." << std::endl;
}
