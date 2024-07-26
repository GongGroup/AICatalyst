import pandas as pd
import numpy as np
from AICatalysis.common.constant import ModelDataDir, ModelData


# 加载数据集
def load_model_dataset(model_dataset=ModelData):
    model_dataset = pd.read_csv(model_dataset, index_col=0)
    return model_dataset


fill_0 = ['cycle', 'c_sub2', 'c_ligand', 'c_base1', 'c_base2', 'c_add', 'c_oxidant', 'pressure CO', 'pressure O2', 'pressure other']
fill_1 = ['solvent_ratio']
fill_mean = ['time']

def fillna(model_dataset):
    model_dataset[fill_0] = model_dataset[fill_0].fillna(0)
    model_dataset[fill_1] = model_dataset[fill_1].fillna(1)
    model_dataset[fill_mean] = model_dataset[fill_mean].fillna(model_dataset[fill_mean].mean())
    model_dataset.iloc[:, 5:] =  model_dataset.iloc[:, 5:].fillna(0)
    return model_dataset


def delta_sub(model_dataset):
    import re
    sub1_col = [col for col in model_dataset.columns if col.startswith('sub1_')]
    sub2_col = [col for col in model_dataset.columns if col.startswith('sub2_')]
    product_col = [col for col in model_dataset.columns if col.startswith('product_')]
    delta_sub1_col = [re.sub('sub1', 'delta_sub1', col) for col in sub1_col]
    delta_sub2_col = [re.sub('sub2', 'delta_sub2', col) for col in sub2_col]
    delta_sub1_pd = pd.DataFrame(model_dataset[product_col].values - model_dataset[sub1_col].values,
                                 columns=delta_sub1_col)
    delta_sub2_pd = pd.DataFrame(model_dataset[product_col].values - model_dataset[sub2_col].values,
                                 columns=delta_sub2_col)
    model_dataset = pd.concat([model_dataset, delta_sub1_pd, delta_sub2_pd], axis=1)
    model_dataset.drop(columns=sub1_col+sub2_col+product_col, inplace=True)
    return model_dataset


def drop_nan_zero(model_dataset):
    nunique_value_cols = model_dataset.columns[model_dataset.nunique() == 1]
    model_dataset.drop(nunique_value_cols, axis=1, inplace=True)
    model_dataset.dropna(how='all', axis=1, inplace=True)
    return model_dataset

def generate_rcol(model_dataset):
    rcol = pd.Series([set([r1, r2])
                      for r1,r2 in zip(model_dataset.r1.values, model_dataset.r2.values)])
    model_dataset.insert(4, 'r', rcol)
    return model_dataset

def move_yield_to_end(model_dataset):
    cols = list(model_dataset.columns)
    cols.append(cols.pop(cols.index('yield')))
    model_dataset = model_dataset.loc[:, cols]
    return model_dataset


def drop_low_std(model_dataset, thresh=0.1):
    dataset_mean = model_dataset.mean()[2:-1]
    dataset_std = model_dataset.std()[2:-1]
    cols_to_drop = dataset_std[dataset_std < thresh * dataset_mean].index
    return model_dataset.drop(cols_to_drop, axis=1)


def standardize(model_dataset):
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    model_dataset.iloc[:, 6: -1] = scaler.fit_transform(model_dataset.iloc[:, 6: -1])
    return model_dataset

def generate_onehot(model_dataset, column):
    onehot_col = pd.get_dummies(model_dataset[column], prefix=column)
    for index, col in enumerate(onehot_col.columns):
        model_dataset.insert(model_dataset.columns.get_loc(column) + index + 1, col, onehot_col[col])
    return model_dataset

def split_model_dataset(model_dataset=ModelData, save=False, drop_reaction_type=False):
    if drop_reaction_type:
        model_dataset_A = model_dataset[model_dataset.reaction_type == 'A'].drop('reaction_type', axis=1)
        model_dataset_B = model_dataset[model_dataset.reaction_type == 'B'].drop('reaction_type', axis=1)
    else:
        model_dataset_A = model_dataset[model_dataset.reaction_type == 'A']
        model_dataset_B = model_dataset[model_dataset.reaction_type == 'B']
    if save:
        model_dataset_A.to_csv(ModelDataDir / 'model_dataset_A.csv', index=False)
        model_dataset_B.to_csv(ModelDataDir / 'model_dataset_B.csv', index=False)
    else:
        return model_dataset_A, model_dataset_B


def model_preprocessing(model_dataset=ModelData, save=False, delta=False):
    model_dataset = load_model_dataset(model_dataset)
    model_dataset = fillna(model_dataset)
    if delta:
        model_dataset = delta_sub(model_dataset)
    model_dataset = drop_nan_zero(model_dataset)
    model_dataset = generate_rcol(model_dataset)
    model_dataset = move_yield_to_end(model_dataset)
    # model_dataset = drop_low_std(model_dataset, thresh=0.1)
    model_dataset = standardize(model_dataset)
    model_dataset = generate_onehot(model_dataset, 'cycle')
    # model_dataset = generate_onehot(model_dataset, 'r')
    if save:
        split_model_dataset(model_dataset, save)
    else:
        model_dataset_A, model_dataset_B = split_model_dataset(model_dataset, save)
        return model_dataset, model_dataset_A, model_dataset_B

# model_preprocessing(save=True, delta=True)

def pca_fit_transform(n_components, X_train, X_test=None):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_components)
    X_train = pca.fit_transform(X_train)
    if X_test is None:
        return X_train
    else:
        X_test = pca.transform(X_test)
        return X_train, X_test

def split_chemical_condition(train_set, test_set, n_components=90):
    from tensorflow import keras
    import tensorflow as tf
    train_condition_data, train_chemical_data = tf.convert_to_tensor(train_set.condition), tf.convert_to_tensor(
        train_set.chemical)
    test_condition_data, test_chemical_data = tf.convert_to_tensor(test_set.condition), tf.convert_to_tensor(
        test_set.chemical)

    X_train_chemical, X_test_chemical = pca_fit_transform(n_components, train_chemical_data, test_chemical_data)

    X_train = np.concatenate((train_condition_data, X_train_chemical), axis=1)
    X_test  = np.concatenate((test_condition_data, X_test_chemical), axis=1)

    y_train = np.array(train_set.y)
    y_test = np.array(test_set.y)

    return X_train, X_test, y_train, y_test

def ss_split(data, label, test_size, random_state, n_splits=1):
    split = StratifiedShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)
    for train_index, test_index in split.split(data, data[label]):
        train_set = data.iloc[train_index]
        test_set = data.iloc[test_index]
    return train_set, test_set

def lasso_gridcv(X_train, y_train, init_model=None):
    from sklearn.linear_model import Lasso
    if init_model == None:
        init_model = Lasso()
    param_grid_list = [{'alpha': [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]}]
    return  grid_search(init_model, param_grid_list, X_train, y_train)

def dt_gridcv(X_train, y_train, init_model=None):
    from sklearn.tree import DecisionTreeRegressor
    if init_model == None:
        init_model = DecisionTreeRegressor()
    param_grid_list = [{'max_depth': list(range(3, 20))},
                       {'min_samples_leaf': list(range(1, 4, 1)), 'max_leaf_nodes': list(range(10, 60, 5))}]
    return grid_search(init_model, param_grid_list, X_train, y_train)


def xgb_gridcv(X_train, y_train, init_model=None):
    import xgboost
    if init_model == None:
        init_model = xgboost.XGBRegressor(eta=0.11, colsample_bytree=0.76, subsample=0.83, max_depth=3, min_child_weight=15,
                                         gamma=0)
    param_grid_list = [{'max_depth': list(range(3,10)), 'min_child_weight': list(range(4,16))},
                       {'gamma': np.arange(0, 0.5, 0.1)},
                       {'eta': np.arange(0.1,0.2,0.01)},
                       {'colsample_bytree': np.arange(0.7, 0.99, 0.03),'subsample': np.arange(0.7,0.99,0.03)}]

    return grid_search(init_model, param_grid_list, X_train, y_train)


def rf_gridcv(X_train, y_train, init_model=None):
    from sklearn.ensemble import RandomForestRegressor
    if init_model == None:
        init_model = RandomForestRegressor(n_estimators=100, max_depth=12, min_samples_leaf=2, max_leaf_nodes=35)
    param_grid_list = [{'max_depth': np.arange(3,20,1)},
                       {'min_samples_leaf': np.arange(1,4,1), 'max_leaf_nodes': np.arange(10,60,5)}]

    return grid_search(init_model, param_grid_list, X_train, y_train)

def svr_gridcv(X_train, y_train, init_model=None):
    from sklearn.svm import SVR
    if init_model == None:
        init_model = SVR(kernel="rbf", gamma=0.01, C=40, epsilon=5.5)
    param_grid_list = [{'gamma': np.arange(0.01,0.5,0.01)},
                       {'C': np.arange(30,100, 10)},
                       {'epsilon': np.arange(2,10, 0.1)}]
    return grid_search(init_model, param_grid_list, X_train, y_train)

def grid_search(init_model, param_grid_list, X_train, y_train):
    from sklearn.model_selection import GridSearchCV
    for i, param_grid in enumerate(param_grid_list):
        if i == 0:
            grid = GridSearchCV(init_model,  param_grid=param_grid, cv=5, scoring='neg_mean_squared_error')
        else:
            grid = GridSearchCV(grid.best_estimator_,
                            param_grid=param_grid, cv=5, scoring='neg_mean_squared_error', verbose=True)
        grid.fit(np.array(X_train), np.array(y_train))
    return grid

