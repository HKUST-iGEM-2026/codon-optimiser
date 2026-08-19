// Deployment shim, not a second implementation.
//
// FastAPI in the container already serves the page, the static files and /api/v1.
// Cloudflare will not route to a container without a Worker in front of it, so this
// Worker proxies every request unchanged and decides nothing. Anything that looks
// like product logic belongs in app/, where the tests can see it.
import { Container } from "@cloudflare/containers";

export class CodonOptimiser extends Container {
  defaultPort = 8000; // the Dockerfile's `uvicorn --port 8000`

  // Container disk is ephemeral: when the instance sleeps, data/permalinks.sqlite
  // is gone and old /p/{id} links 404. sleepAfter defaults to 10 minutes of idle.
  // Raising it only postpones that while paying for an idle instance, so the real
  // fix, if permalinks must outlive an instance, is to move the rows off local disk.

  // The container's own socket peer is this Worker, so without a real client IP the
  // per-IP CPU budget in app/api.py collapses into one shared bucket and a single
  // caller can exhaust everyone's. Cloudflare sets CF-Connecting-IP on the request
  // it forwards; naming the header here is what tells the app to trust it.
  envVars = { CLIENT_IP_HEADER: "CF-Connecting-IP" };

  onError(error) {
    console.error("container error:", error);
  }
}

export default {
  async fetch(request, env) {
    // One named instance, matching max_instances: 1 in wrangler.jsonc.
    return env.OPTIMISER.getByName("app").fetch(request);
  },
};
