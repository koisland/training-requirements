import json
from fastapi import FastAPI, Request, HTTPException
from databases import Database


db = Database("sqlite:///data/gapminder.db")

app = FastAPI()

DB_COLS = {
    "cty": 'Country Name',
    "year": 'Year',
    "agr": 'Agriculture, value added (% of GDP)',
    "co2": 'CO2 emissions (metric tons per capita)',
    "fin": 'Domestic credit provided by financial sector (% of GDP)',
    "elec": 'Electric power consumption (kWh per capita)',
    "energy": 'Energy use (kg of oil equivalent per capita)',
    "exp": 'Exports of goods and services (% of GDP)',
    "fert": 'Fertility rate, total (births per woman)',
    "gdp_g": 'GDP growth (annual %)',
    "imp": 'Imports of goods and services (% of GDP)',
    "ind": 'Industry, value added (% of GDP)',
    "inf": 'Inflation, GDP deflator (annual %)',
    "lexp": 'Life expectancy at birth, total (years)',
    "pden": 'Population density (people per sq. km of land area)',
    "serv": 'Services, etc., value added (% of GDP)',
    "pop": 'pop',
    "cont": 'continent',
    "gdp_c": 'gdpPercap'

}

Q_COND = {
    "gt": ">",
    "ge": ">=",
    "lt": "<",
    "le": "<=",
    "eq": "=",
    "ne": "!="
}


@app.on_event("startup")
async def db_connect():
    await db.connect()


@app.on_event("shutdown")
async def db_shutdown():
    await db.disconnect()


async def construct_sql(params):
    """
    Construct SQL condition strings from given requests.
    @param params: - [dict] of parameters from request.query_params
    @return:
        new_params - [dict] cleaned and converted params
        param_q - [dict] of column names and parameterized SQL condition strings
    """
    # Cleaned params
    new_params = {}
    # Dict column and SQL condition strings
    param_q = {}
    for col, val in params.items():
        # Separate cond from column.
        split_col = col.split("-")
        colname = split_col[0]

        # Convert values to numeric if possible. otherwise, capitalize and replace hyphen with space
        if num_val := val.isnumeric():
            new_params[colname] = float(val)
        else:
            new_params[colname] = val.title().replace("-", " ")

        if colname in DB_COLS:
            if len(split_col) == 2:
                # Construct SQL condition. colname and cond are safe since must be in above dicts. val sanitized below.
                # Get the condition if in listed conditions. Otherwise, =.
                # TODO: Check if condtion make sense for col.
                param_q[colname] = f'[{DB_COLS[colname]}] {Q_COND.get(split_col[1], "=")} :{colname}'
            else:
                param_q[colname] = f'[{DB_COLS[colname]}] = :{colname}'
        else:
            raise HTTPException(status_code=404, detail="Column not found")

    return param_q, new_params


@app.get("/api/country")
async def country(request: Request):
    params = request.query_params
    param_q, new_params = await construct_sql(params)

    query = f"SELECT [Country Name] FROM gapminder WHERE {' AND '.join(param_q.values())}"
    res = await db.fetch_all(query=query, values=new_params)

    return {"countries": res}


@app.get("/api/gapminder")
async def gapminder(request: Request):
    params = request.query_params
    param_q, new_params = await construct_sql(params)

    query = f"SELECT * FROM gapminder WHERE {' AND '.join(param_q.values())}"
    res = await db.fetch_all(query=query, values=new_params)

    # Convert row object to list of dict and dump as json.
    return json.dumps([dict(row) for row in res])


# To run:
# uvicorn main:app --reload
